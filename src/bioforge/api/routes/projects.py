"""Projects routes."""

from fastapi import APIRouter, HTTPException, status
from pydantic import BaseModel
from sqlalchemy import select

from bioforge.api.deps import SessionDep, CurrentActiveUser
from bioforge.models import Project

router = APIRouter(prefix="/projects", tags=["projects"])


# ─────────────────────────────────────────────────────────────────────────────
# Schemas
# ─────────────────────────────────────────────────────────────────────────────

class ProjectCreate(BaseModel):
    name: str
    description: str | None = None
    reference_genome: str = "GRCh38"


class ProjectUpdate(BaseModel):
    name: str | None = None
    description: str | None = None
    reference_genome: str | None = None


class ProjectResponse(BaseModel):
    id: int
    name: str
    description: str | None
    reference_genome: str
    owner_id: int

    model_config = {"from_attributes": True}


# ─────────────────────────────────────────────────────────────────────────────
# Routes
# ─────────────────────────────────────────────────────────────────────────────

@router.get("", response_model=list[ProjectResponse])
async def list_projects(session: SessionDep, user: CurrentActiveUser):
    """List all projects for current user."""
    result = await session.execute(
        select(Project).where(Project.owner_id == user.id)
    )
    return result.scalars().all()


@router.post("", response_model=ProjectResponse, status_code=status.HTTP_201_CREATED)
async def create_project(
    data: ProjectCreate,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Create a new project."""
    project = Project(
        name=data.name,
        description=data.description,
        reference_genome=data.reference_genome,
        owner_id=user.id,
    )
    session.add(project)
    await session.commit()
    await session.refresh(project)
    return project


@router.get("/{project_id}", response_model=ProjectResponse)
async def get_project(
    project_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Get project by ID."""
    result = await session.execute(
        select(Project).where(
            Project.id == project_id,
            Project.owner_id == user.id,
        )
    )
    project = result.scalar_one_or_none()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    return project


@router.patch("/{project_id}", response_model=ProjectResponse)
async def update_project(
    project_id: int,
    data: ProjectUpdate,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Update project."""
    result = await session.execute(
        select(Project).where(
            Project.id == project_id,
            Project.owner_id == user.id,
        )
    )
    project = result.scalar_one_or_none()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    
    update_data = data.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(project, field, value)
    
    await session.commit()
    await session.refresh(project)
    return project


@router.delete("/{project_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_project(
    project_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Delete project."""
    result = await session.execute(
        select(Project).where(
            Project.id == project_id,
            Project.owner_id == user.id,
        )
    )
    project = result.scalar_one_or_none()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    
    await session.delete(project)
    await session.commit()
