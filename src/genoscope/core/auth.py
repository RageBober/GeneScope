"""
Authentication and Authorization Module for GenoScope

Provides JWT tokens, OAuth2, API keys, and role-based access control.
"""

from __future__ import annotations

import hashlib
import hmac
import logging
import secrets
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Tuple
from functools import wraps

import jwt
from passlib.context import CryptContext
from pydantic import BaseModel, EmailStr, Field, validator

logger = logging.getLogger(__name__)

# Password hashing context
pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

# JWT settings
JWT_SECRET_KEY = "your-secret-key-change-in-production"  # Should be from env
JWT_ALGORITHM = "HS256"
JWT_ACCESS_TOKEN_EXPIRE_MINUTES = 30
JWT_REFRESH_TOKEN_EXPIRE_DAYS = 7


class UserRole(Enum):
    """User roles with hierarchical permissions."""
    
    ADMIN = "admin"  # Full system access
    RESEARCHER = "researcher"  # Can run analyses, view all data
    COLLABORATOR = "collaborator"  # Can view shared data
    VIEWER = "viewer"  # Read-only access
    API_USER = "api_user"  # Programmatic access


class Permission(Enum):
    """System permissions."""
    
    # Analysis permissions
    CREATE_ANALYSIS = "create_analysis"
    VIEW_ANALYSIS = "view_analysis"
    DELETE_ANALYSIS = "delete_analysis"
    SHARE_ANALYSIS = "share_analysis"
    
    # Data permissions
    UPLOAD_DATA = "upload_data"
    VIEW_DATA = "view_data"
    DELETE_DATA = "delete_data"
    EXPORT_DATA = "export_data"
    
    # User management
    CREATE_USER = "create_user"
    MODIFY_USER = "modify_user"
    DELETE_USER = "delete_user"
    VIEW_USERS = "view_users"
    
    # System permissions
    VIEW_LOGS = "view_logs"
    MODIFY_SETTINGS = "modify_settings"
    VIEW_BILLING = "view_billing"
    MANAGE_BILLING = "manage_billing"


# Role-permission mapping
ROLE_PERMISSIONS: Dict[UserRole, Set[Permission]] = {
    UserRole.ADMIN: set(Permission),  # All permissions
    UserRole.RESEARCHER: {
        Permission.CREATE_ANALYSIS,
        Permission.VIEW_ANALYSIS,
        Permission.DELETE_ANALYSIS,
        Permission.SHARE_ANALYSIS,
        Permission.UPLOAD_DATA,
        Permission.VIEW_DATA,
        Permission.DELETE_DATA,
        Permission.EXPORT_DATA,
        Permission.VIEW_BILLING,
    },
    UserRole.COLLABORATOR: {
        Permission.VIEW_ANALYSIS,
        Permission.VIEW_DATA,
        Permission.EXPORT_DATA,
    },
    UserRole.VIEWER: {
        Permission.VIEW_ANALYSIS,
        Permission.VIEW_DATA,
    },
    UserRole.API_USER: {
        Permission.CREATE_ANALYSIS,
        Permission.VIEW_ANALYSIS,
        Permission.UPLOAD_DATA,
        Permission.VIEW_DATA,
        Permission.EXPORT_DATA,
    },
}


@dataclass
class User:
    """User model."""
    
    id: str
    email: str
    username: str
    full_name: str
    role: UserRole
    is_active: bool = True
    is_verified: bool = False
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    last_login: Optional[datetime] = None
    password_hash: Optional[str] = None
    api_keys: List[str] = field(default_factory=list)
    oauth_providers: Dict[str, str] = field(default_factory=dict)  # provider -> user_id
    permissions: Set[Permission] = field(default_factory=set)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Initialize permissions based on role."""
        if not self.permissions:
            self.permissions = ROLE_PERMISSIONS.get(self.role, set())
    
    def has_permission(self, permission: Permission) -> bool:
        """Check if user has a specific permission."""
        return permission in self.permissions
    
    def can_access_resource(self, resource_owner_id: str) -> bool:
        """Check if user can access a resource."""
        # Admins can access everything
        if self.role == UserRole.ADMIN:
            return True
        
        # Users can access their own resources
        if resource_owner_id == self.id:
            return True
        
        # Check if resource is shared with user (implement sharing logic)
        # This would check a sharing database/table
        return False


class TokenData(BaseModel):
    """JWT token payload model."""
    
    sub: str  # Subject (user ID)
    email: Optional[str] = None
    role: Optional[str] = None
    permissions: List[str] = []
    exp: Optional[datetime] = None
    iat: Optional[datetime] = None
    jti: Optional[str] = None  # JWT ID for revocation


class LoginRequest(BaseModel):
    """Login request model."""
    
    email: EmailStr
    password: str


class RegisterRequest(BaseModel):
    """Registration request model."""
    
    email: EmailStr
    username: str = Field(..., min_length=3, max_length=50)
    password: str = Field(..., min_length=8)
    full_name: str
    
    @validator("password")
    def validate_password(cls, v):
        """Ensure password meets security requirements."""
        if not any(c.isupper() for c in v):
            raise ValueError("Password must contain at least one uppercase letter")
        if not any(c.islower() for c in v):
            raise ValueError("Password must contain at least one lowercase letter")
        if not any(c.isdigit() for c in v):
            raise ValueError("Password must contain at least one digit")
        return v


class AuthManager:
    """Manages authentication and authorization."""
    
    def __init__(self, secret_key: Optional[str] = None):
        """
        Initialize auth manager.
        
        Args:
            secret_key: JWT secret key (uses default if not provided)
        """
        self.secret_key = secret_key or JWT_SECRET_KEY
        self.revoked_tokens: Set[str] = set()  # In production, use Redis
        self.api_keys: Dict[str, str] = {}  # api_key -> user_id
        self.users: Dict[str, User] = {}  # user_id -> User (in production, use database)
        
        # OAuth2 providers configuration
        self.oauth_providers = {
            "google": {
                "client_id": "",  # From environment
                "client_secret": "",
                "redirect_uri": "http://localhost:8000/auth/google/callback",
                "auth_url": "https://accounts.google.com/o/oauth2/v2/auth",
                "token_url": "https://oauth2.googleapis.com/token",
                "user_info_url": "https://www.googleapis.com/oauth2/v2/userinfo",
            },
            "github": {
                "client_id": "",
                "client_secret": "",
                "redirect_uri": "http://localhost:8000/auth/github/callback",
                "auth_url": "https://github.com/login/oauth/authorize",
                "token_url": "https://github.com/login/oauth/access_token",
                "user_info_url": "https://api.github.com/user",
            },
        }
    
    def hash_password(self, password: str) -> str:
        """Hash a password."""
        return pwd_context.hash(password)
    
    def verify_password(self, plain_password: str, hashed_password: str) -> bool:
        """Verify a password against its hash."""
        return pwd_context.verify(plain_password, hashed_password)
    
    def create_access_token(
        self,
        user: User,
        expires_delta: Optional[timedelta] = None
    ) -> str:
        """
        Create a JWT access token.
        
        Args:
            user: User object
            expires_delta: Token expiration time
            
        Returns:
            JWT token string
        """
        if expires_delta:
            expire = datetime.now(timezone.utc) + expires_delta
        else:
            expire = datetime.now(timezone.utc) + timedelta(
                minutes=JWT_ACCESS_TOKEN_EXPIRE_MINUTES
            )
        
        # Create token payload
        payload = {
            "sub": user.id,
            "email": user.email,
            "role": user.role.value,
            "permissions": [p.value for p in user.permissions],
            "exp": expire,
            "iat": datetime.now(timezone.utc),
            "jti": secrets.token_urlsafe(16),  # Unique token ID
        }
        
        # Encode token
        token = jwt.encode(payload, self.secret_key, algorithm=JWT_ALGORITHM)
        return token
    
    def create_refresh_token(self, user: User) -> str:
        """
        Create a JWT refresh token.
        
        Args:
            user: User object
            
        Returns:
            JWT refresh token string
        """
        expire = datetime.now(timezone.utc) + timedelta(
            days=JWT_REFRESH_TOKEN_EXPIRE_DAYS
        )
        
        payload = {
            "sub": user.id,
            "type": "refresh",
            "exp": expire,
            "iat": datetime.now(timezone.utc),
            "jti": secrets.token_urlsafe(16),
        }
        
        token = jwt.encode(payload, self.secret_key, algorithm=JWT_ALGORITHM)
        return token
    
    def verify_token(self, token: str) -> Optional[TokenData]:
        """
        Verify and decode a JWT token.
        
        Args:
            token: JWT token string
            
        Returns:
            TokenData if valid, None otherwise
        """
        try:
            # Decode token
            payload = jwt.decode(
                token,
                self.secret_key,
                algorithms=[JWT_ALGORITHM]
            )
            
            # Check if token is revoked
            jti = payload.get("jti")
            if jti and jti in self.revoked_tokens:
                logger.warning(f"Attempted to use revoked token: {jti}")
                return None
            
            # Create TokenData
            return TokenData(**payload)
            
        except jwt.ExpiredSignatureError:
            logger.debug("Token has expired")
            return None
        except jwt.InvalidTokenError as e:
            logger.debug(f"Invalid token: {e}")
            return None
    
    def revoke_token(self, token: str) -> bool:
        """
        Revoke a JWT token.
        
        Args:
            token: JWT token to revoke
            
        Returns:
            True if revoked successfully
        """
        try:
            payload = jwt.decode(
                token,
                self.secret_key,
                algorithms=[JWT_ALGORITHM],
                options={"verify_exp": False}  # Allow expired tokens to be revoked
            )
            
            jti = payload.get("jti")
            if jti:
                self.revoked_tokens.add(jti)
                logger.info(f"Token revoked: {jti}")
                return True
                
        except jwt.InvalidTokenError:
            pass
        
        return False
    
    def create_api_key(self, user: User) -> str:
        """
        Create an API key for a user.
        
        Args:
            user: User object
            
        Returns:
            API key string
        """
        # Generate secure random key
        api_key = f"gsk_{secrets.token_urlsafe(32)}"  # gsk = genoscope key
        
        # Store mapping
        self.api_keys[api_key] = user.id
        user.api_keys.append(api_key)
        
        logger.info(f"Created API key for user {user.id}")
        return api_key
    
    def verify_api_key(self, api_key: str) -> Optional[User]:
        """
        Verify an API key.
        
        Args:
            api_key: API key to verify
            
        Returns:
            User if valid, None otherwise
        """
        user_id = self.api_keys.get(api_key)
        if user_id:
            user = self.users.get(user_id)
            if user and user.is_active:
                return user
        
        return None
    
    def revoke_api_key(self, api_key: str) -> bool:
        """
        Revoke an API key.
        
        Args:
            api_key: API key to revoke
            
        Returns:
            True if revoked successfully
        """
        if api_key in self.api_keys:
            user_id = self.api_keys[api_key]
            del self.api_keys[api_key]
            
            # Remove from user's api_keys list
            user = self.users.get(user_id)
            if user and api_key in user.api_keys:
                user.api_keys.remove(api_key)
            
            logger.info(f"API key revoked: {api_key[:10]}...")
            return True
        
        return False
    
    def register_user(self, request: RegisterRequest) -> Tuple[User, str]:
        """
        Register a new user.
        
        Args:
            request: Registration request
            
        Returns:
            Tuple of (User, access_token)
        """
        # Check if user exists
        for user in self.users.values():
            if user.email == request.email:
                raise ValueError("Email already registered")
            if user.username == request.username:
                raise ValueError("Username already taken")
        
        # Create new user
        user = User(
            id=secrets.token_urlsafe(16),
            email=request.email,
            username=request.username,
            full_name=request.full_name,
            role=UserRole.RESEARCHER,  # Default role
            password_hash=self.hash_password(request.password),
        )
        
        # Store user
        self.users[user.id] = user
        
        # Create access token
        access_token = self.create_access_token(user)
        
        logger.info(f"Registered new user: {user.email}")
        return user, access_token
    
    def login(self, request: LoginRequest) -> Tuple[User, str, str]:
        """
        Login a user.
        
        Args:
            request: Login request
            
        Returns:
            Tuple of (User, access_token, refresh_token)
        """
        # Find user by email
        user = None
        for u in self.users.values():
            if u.email == request.email:
                user = u
                break
        
        if not user:
            raise ValueError("Invalid email or password")
        
        # Verify password
        if not user.password_hash or not self.verify_password(
            request.password, user.password_hash
        ):
            raise ValueError("Invalid email or password")
        
        # Check if user is active
        if not user.is_active:
            raise ValueError("User account is disabled")
        
        # Update last login
        user.last_login = datetime.now(timezone.utc)
        
        # Create tokens
        access_token = self.create_access_token(user)
        refresh_token = self.create_refresh_token(user)
        
        logger.info(f"User logged in: {user.email}")
        return user, access_token, refresh_token
    
    def refresh_access_token(self, refresh_token: str) -> Optional[str]:
        """
        Refresh an access token using a refresh token.
        
        Args:
            refresh_token: Refresh token
            
        Returns:
            New access token if valid, None otherwise
        """
        token_data = self.verify_token(refresh_token)
        if not token_data:
            return None
        
        # Check if it's a refresh token
        if token_data.role:  # Access tokens have role
            return None
        
        # Get user
        user = self.users.get(token_data.sub)
        if not user or not user.is_active:
            return None
        
        # Create new access token
        return self.create_access_token(user)
    
    def get_oauth_redirect_url(self, provider: str, state: str) -> str:
        """
        Get OAuth2 redirect URL.
        
        Args:
            provider: OAuth provider (google, github)
            state: State parameter for CSRF protection
            
        Returns:
            OAuth redirect URL
        """
        if provider not in self.oauth_providers:
            raise ValueError(f"Unknown OAuth provider: {provider}")
        
        config = self.oauth_providers[provider]
        
        params = {
            "client_id": config["client_id"],
            "redirect_uri": config["redirect_uri"],
            "response_type": "code",
            "state": state,
        }
        
        if provider == "google":
            params["scope"] = "openid email profile"
        elif provider == "github":
            params["scope"] = "user:email"
        
        # Build URL
        from urllib.parse import urlencode
        return f"{config['auth_url']}?{urlencode(params)}"
    
    def handle_oauth_callback(
        self,
        provider: str,
        code: str,
        state: str
    ) -> Tuple[User, str, str]:
        """
        Handle OAuth2 callback.
        
        Args:
            provider: OAuth provider
            code: Authorization code
            state: State parameter
            
        Returns:
            Tuple of (User, access_token, refresh_token)
        """
        # Exchange code for token
        # Get user info from provider
        # Create or update user
        # Return tokens
        
        # This is a simplified implementation
        # In production, make actual HTTP requests to provider
        raise NotImplementedError("OAuth callback handling not fully implemented")


# Decorators for FastAPI/Flask
def require_auth(permission: Optional[Permission] = None):
    """
    Decorator to require authentication and optionally check permission.
    
    Args:
        permission: Required permission (optional)
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get token from request (implementation depends on framework)
            # Verify token
            # Check permission if specified
            # Call function with user context
            return func(*args, **kwargs)
        return wrapper
    return decorator


# Session manager for web sessions
class SessionManager:
    """Manages web sessions with Redis backend."""
    
    def __init__(self, redis_client=None):
        """
        Initialize session manager.
        
        Args:
            redis_client: Redis client instance
        """
        self.redis = redis_client
        self.sessions: Dict[str, Dict[str, Any]] = {}  # Fallback to memory if no Redis
    
    def create_session(self, user_id: str, data: Dict[str, Any]) -> str:
        """Create a new session."""
        session_id = secrets.token_urlsafe(32)
        session_data = {
            "user_id": user_id,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "last_activity": datetime.now(timezone.utc).isoformat(),
            **data
        }
        
        if self.redis:
            # Store in Redis with TTL
            self.redis.setex(
                f"session:{session_id}",
                timedelta(hours=24),
                json.dumps(session_data)
            )
        else:
            # Store in memory
            self.sessions[session_id] = session_data
        
        return session_id
    
    def get_session(self, session_id: str) -> Optional[Dict[str, Any]]:
        """Get session data."""
        if self.redis:
            data = self.redis.get(f"session:{session_id}")
            if data:
                return json.loads(data)
        else:
            return self.sessions.get(session_id)
        
        return None
    
    def update_session(self, session_id: str, data: Dict[str, Any]) -> bool:
        """Update session data."""
        session = self.get_session(session_id)
        if not session:
            return False
        
        session.update(data)
        session["last_activity"] = datetime.now(timezone.utc).isoformat()
        
        if self.redis:
            self.redis.setex(
                f"session:{session_id}",
                timedelta(hours=24),
                json.dumps(session)
            )
        else:
            self.sessions[session_id] = session
        
        return True
    
    def delete_session(self, session_id: str) -> bool:
        """Delete a session."""
        if self.redis:
            return bool(self.redis.delete(f"session:{session_id}"))
        else:
            if session_id in self.sessions:
                del self.sessions[session_id]
                return True
        
        return False


# Rate limiter
class RateLimiter:
    """Rate limiter for API endpoints."""
    
    def __init__(self, redis_client=None):
        """
        Initialize rate limiter.
        
        Args:
            redis_client: Redis client instance
        """
        self.redis = redis_client
        self.limits: Dict[str, List[float]] = {}  # Fallback to memory
    
    def is_allowed(
        self,
        key: str,
        limit: int = 100,
        window: int = 3600
    ) -> bool:
        """
        Check if request is allowed under rate limit.
        
        Args:
            key: Rate limit key (e.g., user_id, ip_address)
            limit: Maximum requests in window
            window: Time window in seconds
            
        Returns:
            True if allowed, False if rate limited
        """
        current_time = time.time()
        
        if self.redis:
            # Use Redis for distributed rate limiting
            pipe = self.redis.pipeline()
            redis_key = f"rate_limit:{key}"
            
            # Remove old entries
            pipe.zremrangebyscore(redis_key, 0, current_time - window)
            
            # Count current entries
            pipe.zcard(redis_key)
            
            # Add current request
            pipe.zadd(redis_key, {str(current_time): current_time})
            
            # Set expiry
            pipe.expire(redis_key, window)
            
            results = pipe.execute()
            count = results[1]
            
            return count < limit
        else:
            # In-memory rate limiting
            if key not in self.limits:
                self.limits[key] = []
            
            # Remove old entries
            self.limits[key] = [
                t for t in self.limits[key]
                if t > current_time - window
            ]
            
            # Check limit
            if len(self.limits[key]) >= limit:
                return False
            
            # Add current request
            self.limits[key].append(current_time)
            return True


# Audit logger
class AuditLogger:
    """Logs all authentication and authorization events."""
    
    def __init__(self, logger_name: str = "audit"):
        """Initialize audit logger."""
        self.logger = logging.getLogger(logger_name)
    
    def log_login(self, user_id: str, ip_address: str, success: bool):
        """Log login attempt."""
        self.logger.info(
            f"LOGIN {'SUCCESS' if success else 'FAILED'} - "
            f"User: {user_id}, IP: {ip_address}"
        )
    
    def log_access(
        self,
        user_id: str,
        resource: str,
        action: str,
        allowed: bool
    ):
        """Log resource access attempt."""
        self.logger.info(
            f"ACCESS {'GRANTED' if allowed else 'DENIED'} - "
            f"User: {user_id}, Resource: {resource}, Action: {action}"
        )
    
    def log_api_call(
        self,
        user_id: str,
        endpoint: str,
        method: str,
        status_code: int
    ):
        """Log API call."""
        self.logger.info(
            f"API CALL - User: {user_id}, "
            f"Endpoint: {endpoint}, Method: {method}, "
            f"Status: {status_code}"
        )
