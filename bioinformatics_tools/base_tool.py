"""
🔧 Base Class for Bioinformatics Tools

Provides:
- Unified interface for all tools
- Docker containerization
- Input validation
- Resource management
- Error handling
- Logging
"""

from __future__ import annotations

import json
import logging
import shlex
import subprocess
import tempfile
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional

from .config import ToolConfig
from .validators import (
    ValidationError,
    validate_command_args,
    validate_docker_image_name,
    validate_file_path,
    validate_output_directory,
)

logger = logging.getLogger(__name__)


class ToolExecutionError(Exception):
    """Raised when tool execution fails"""
    pass


class BaseBioTool(ABC):
    """
    Abstract base class for bioinformatics tools.

    All tools must implement:
    - _build_command(): Construct command-line arguments
    - _parse_output(): Parse tool output into structured format

    Security features:
    - Input validation
    - Docker isolation
    - Resource limiting
    - Command injection prevention
    """

    def __init__(self, config: ToolConfig):
        """
        Initialize tool with configuration.

        Args:
            config: Tool configuration object
        """
        self.config = config
        self.logger = logging.getLogger(f"{__name__}.{config.name}")

    @abstractmethod
    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """
        Build command-line arguments for the tool.

        Args:
            input_path: Validated input file path
            output_path: Validated output path
            **kwargs: Additional tool-specific arguments

        Returns:
            List of command arguments
        """
        pass

    @abstractmethod
    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """
        Parse tool output into structured format.

        Args:
            output_path: Path to tool output file

        Returns:
            Dictionary with parsed results
        """
        pass

    def validate_inputs(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
    ) -> tuple[Path, Path]:
        """
        Validate input file and output directory.

        Args:
            input_path: Input file path
            output_dir: Output directory (default: config.work_dir)

        Returns:
            Tuple of (validated_input_path, validated_output_dir)

        Raises:
            ValidationError: If validation fails
        """
        # Validate input file
        input_file = validate_file_path(input_path, must_exist=True)
        self.logger.info(f"Validated input: {input_file}")

        # Validate output directory
        if output_dir is None:
            output_dir = self.config.work_dir

        output_directory = validate_output_directory(
            output_dir,
            create_if_missing=True,
        )
        self.logger.info(f"Validated output directory: {output_directory}")

        return input_file, output_directory

    def _build_docker_command(
        self,
        tool_command: list[str],
        input_path: Path,
        output_dir: Path,
    ) -> list[str]:
        """
        Wrap tool command in Docker container.

        Security features:
        - Read-only root filesystem
        - No new privileges
        - Non-root user
        - CPU/memory limits
        - No network access

        Args:
            tool_command: Tool command to run
            input_path: Input file path
            output_dir: Output directory path

        Returns:
            Docker command list
        """
        docker_cfg = self.config.docker_config

        # Validate Docker image name
        if not self.config.docker_image:
            raise ValueError(f"No Docker image specified for {self.config.name}")

        image = validate_docker_image_name(self.config.docker_image)

        # Build Docker command
        cmd = [
            "docker", "run",
            "--rm",  # Remove container after exit
            "--user", docker_cfg.user,  # Non-root user
            "--network", docker_cfg.network_mode,  # No network by default
        ]

        # Resource limits
        cmd.extend([
            "--cpus", str(docker_cfg.cpu_count),
            "--cpu-quota", str(docker_cfg.cpu_quota),
            "--memory", docker_cfg.memory_limit,
            "--memory-swap", docker_cfg.memory_swap_limit,
        ])

        # Security options
        if docker_cfg.read_only_root:
            cmd.append("--read-only")

        if docker_cfg.no_new_privileges:
            cmd.append("--security-opt=no-new-privileges")

        # Mount volumes (read-only input, writable output)
        input_mount = f"{input_path.parent}:/input:ro"
        output_mount = f"{output_dir}:/output:rw"

        cmd.extend(["-v", input_mount, "-v", output_mount])

        # Working directory
        cmd.extend(["-w", "/output"])

        # Docker image
        cmd.append(image)

        # Tool command
        cmd.extend(tool_command)

        return cmd

    def execute(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Execute bioinformatics tool.

        Args:
            input_path: Input file path
            output_dir: Output directory
            **kwargs: Tool-specific arguments

        Returns:
            Dictionary with execution results

        Raises:
            ValidationError: If input validation fails
            ToolExecutionError: If tool execution fails
        """
        # Step 1: Validate inputs
        try:
            input_file, output_directory = self.validate_inputs(input_path, output_dir)
        except ValidationError as exc:
            self.logger.error(f"Input validation failed: {exc}")
            raise

        # Step 2: Generate output filename
        output_file = output_directory / f"{input_file.stem}_{self.config.name}_output.txt"

        # Step 3: Build tool command
        try:
            tool_command = self._build_command(input_file, output_file, **kwargs)
            tool_command = validate_command_args(tool_command)
        except Exception as exc:
            self.logger.error(f"Failed to build command: {exc}")
            raise ToolExecutionError(f"Command build failed: {exc}") from exc

        # Step 4: Wrap in Docker if enabled
        if self.config.use_docker:
            command = self._build_docker_command(tool_command, input_file, output_directory)
        else:
            command = tool_command

        # Step 5: Dry run mode (don't execute)
        if self.config.dry_run:
            self.logger.info(f"[DRY RUN] Would execute: {' '.join(command)}")
            return {
                "status": "dry_run",
                "command": command,
                "output_file": str(output_file),
            }

        # Step 6: Execute command
        self.logger.info(f"Executing: {self.config.name}")
        self.logger.debug(f"Command: {' '.join(command)}")

        start_time = time.time()

        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=self.config.docker_config.timeout_seconds,
                check=True,
            )

            execution_time = time.time() - start_time
            self.logger.info(f"Execution completed in {execution_time:.2f}s")

            # Step 7: Parse output
            parsed_output = self._parse_output(output_file)

            return {
                "status": "success",
                "execution_time": execution_time,
                "output_file": str(output_file),
                "results": parsed_output,
                "stdout": result.stdout,
                "stderr": result.stderr,
            }

        except subprocess.TimeoutExpired as exc:
            self.logger.error(f"Execution timeout after {exc.timeout}s")
            raise ToolExecutionError(f"Timeout after {exc.timeout}s") from exc

        except subprocess.CalledProcessError as exc:
            self.logger.error(f"Execution failed with exit code {exc.returncode}")
            self.logger.error(f"stderr: {exc.stderr}")
            raise ToolExecutionError(
                f"Tool failed with exit code {exc.returncode}: {exc.stderr}"
            ) from exc

        except Exception as exc:
            self.logger.exception(f"Unexpected error: {exc}")
            raise ToolExecutionError(f"Unexpected error: {exc}") from exc

    def execute_async(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> subprocess.Popen:
        """
        Execute tool asynchronously (non-blocking).

        Args:
            input_path: Input file path
            output_dir: Output directory
            **kwargs: Tool-specific arguments

        Returns:
            Popen process object
        """
        input_file, output_directory = self.validate_inputs(input_path, output_dir)
        output_file = output_directory / f"{input_file.stem}_{self.config.name}_output.txt"

        tool_command = self._build_command(input_file, output_file, **kwargs)
        tool_command = validate_command_args(tool_command)

        if self.config.use_docker:
            command = self._build_docker_command(tool_command, input_file, output_directory)
        else:
            command = tool_command

        self.logger.info(f"Starting async execution: {self.config.name}")

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        return process

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.config.name}, version={self.config.version})"
