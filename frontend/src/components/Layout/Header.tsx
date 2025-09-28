import React, { useState, useContext } from 'react';
import {
  AppBar,
  Toolbar,
  IconButton,
  Typography,
  Badge,
  Avatar,
  Menu,
  MenuItem,
  Box,
  Tooltip,
  Divider,
  ListItemIcon,
  ListItemText,
  Chip,
  TextField,
  InputAdornment,
} from '@mui/material';
import { useNavigate } from 'react-router-dom';
import {
  Menu as MenuIcon,
  Bell,
  Search,
  Settings,
  LogOut,
  User,
  HelpCircle,
  Activity,
  FileText,
  Download,
} from 'lucide-react';
import { UserContext } from '../../App';

interface HeaderProps {
  onMenuClick: () => void;
  onLogout: () => void;
}

const Header: React.FC<HeaderProps> = ({ onMenuClick, onLogout }) => {
  const navigate = useNavigate();
  const { user } = useContext(UserContext);
  const [anchorElUser, setAnchorElUser] = useState<null | HTMLElement>(null);
  const [anchorElNotifications, setAnchorElNotifications] = useState<null | HTMLElement>(null);
  const [searchOpen, setSearchOpen] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');

  const handleOpenUserMenu = (event: React.MouseEvent<HTMLElement>) => {
    setAnchorElUser(event.currentTarget);
  };

  const handleCloseUserMenu = () => {
    setAnchorElUser(null);
  };

  const handleOpenNotifications = (event: React.MouseEvent<HTMLElement>) => {
    setAnchorElNotifications(event.currentTarget);
  };

  const handleCloseNotifications = () => {
    setAnchorElNotifications(null);
  };

  const handleSearch = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && searchQuery.trim()) {
      navigate(`/search?q=${encodeURIComponent(searchQuery)}`);
      setSearchQuery('');
      setSearchOpen(false);
    }
  };

  const notifications = [
    {
      id: 1,
      title: 'Analysis Complete',
      message: 'RNA-Seq analysis for Sample_001 has completed',
      time: '5 minutes ago',
      type: 'success',
    },
    {
      id: 2,
      title: 'New Report Available',
      message: 'Quality control report ready for Project Alpha',
      time: '1 hour ago',
      type: 'info',
    },
    {
      id: 3,
      title: 'Storage Warning',
      message: 'You have used 85% of your storage quota',
      time: '2 hours ago',
      type: 'warning',
    },
  ];

  return (
    <AppBar
      position="fixed"
      sx={{
        zIndex: (theme) => theme.zIndex.drawer + 1,
        backgroundColor: 'white',
        color: 'text.primary',
        boxShadow: '0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24)',
      }}
    >
      <Toolbar>
        <IconButton
          color="inherit"
          aria-label="open drawer"
          edge="start"
          onClick={onMenuClick}
          sx={{ mr: 2 }}
        >
          <MenuIcon />
        </IconButton>

        <Box sx={{ display: 'flex', alignItems: 'center', flexGrow: 1 }}>
          <Typography variant="h6" noWrap component="div" sx={{ fontWeight: 700 }}>
            GenoScope
          </Typography>
          <Chip
            label="BETA"
            size="small"
            color="primary"
            sx={{ ml: 1 }}
          />
        </Box>

        {/* Search Bar */}
        <Box sx={{ flexGrow: 1, mx: 4, maxWidth: 500 }}>
          <TextField
            size="small"
            fullWidth
            placeholder="Search analyses, samples, or projects..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            onKeyPress={handleSearch}
            InputProps={{
              startAdornment: (
                <InputAdornment position="start">
                  <Search size={18} />
                </InputAdornment>
              ),
              sx: {
                borderRadius: 2,
                backgroundColor: 'grey.50',
              },
            }}
          />
        </Box>

        {/* Activity Monitor */}
        <Tooltip title="System Activity">
          <Box
            sx={{
              display: 'flex',
              alignItems: 'center',
              mr: 2,
              px: 2,
              py: 0.5,
              borderRadius: 1,
              backgroundColor: 'success.light',
              color: 'success.main',
            }}
          >
            <Activity size={16} />
            <Typography variant="caption" sx={{ ml: 1, fontWeight: 600 }}>
              3 Active Jobs
            </Typography>
          </Box>
        </Tooltip>

        {/* Notifications */}
        <IconButton
          size="large"
          color="inherit"
          onClick={handleOpenNotifications}
          sx={{ mr: 1 }}
        >
          <Badge badgeContent={notifications.length} color="error">
            <Bell />
          </Badge>
        </IconButton>

        {/* User Menu */}
        <Box sx={{ flexGrow: 0 }}>
          <Tooltip title="Account settings">
            <IconButton onClick={handleOpenUserMenu} sx={{ p: 0 }}>
              <Avatar
                alt={user?.name || 'User'}
                sx={{ bgcolor: 'primary.main' }}
              >
                {user?.name?.[0] || 'U'}
              </Avatar>
            </IconButton>
          </Tooltip>
        </Box>

        {/* Notifications Menu */}
        <Menu
          sx={{ mt: '45px', maxWidth: 400 }}
          id="notifications-menu"
          anchorEl={anchorElNotifications}
          anchorOrigin={{
            vertical: 'top',
            horizontal: 'right',
          }}
          keepMounted
          transformOrigin={{
            vertical: 'top',
            horizontal: 'right',
          }}
          open={Boolean(anchorElNotifications)}
          onClose={handleCloseNotifications}
        >
          <Box sx={{ px: 2, py: 1 }}>
            <Typography variant="h6">Notifications</Typography>
          </Box>
          <Divider />
          {notifications.map((notification) => (
            <MenuItem
              key={notification.id}
              onClick={handleCloseNotifications}
              sx={{ py: 2, px: 2 }}
            >
              <Box>
                <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>
                  {notification.title}
                </Typography>
                <Typography variant="body2" color="text.secondary">
                  {notification.message}
                </Typography>
                <Typography variant="caption" color="text.secondary">
                  {notification.time}
                </Typography>
              </Box>
            </MenuItem>
          ))}
          <Divider />
          <MenuItem onClick={() => navigate('/notifications')}>
            <Typography variant="body2" color="primary" sx={{ mx: 'auto' }}>
              View All Notifications
            </Typography>
          </MenuItem>
        </Menu>

        {/* User Menu */}
        <Menu
          sx={{ mt: '45px' }}
          id="user-menu"
          anchorEl={anchorElUser}
          anchorOrigin={{
            vertical: 'top',
            horizontal: 'right',
          }}
          keepMounted
          transformOrigin={{
            vertical: 'top',
            horizontal: 'right',
          }}
          open={Boolean(anchorElUser)}
          onClose={handleCloseUserMenu}
        >
          <Box sx={{ px: 2, py: 1 }}>
            <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
              {user?.name || 'Guest User'}
            </Typography>
            <Typography variant="body2" color="text.secondary">
              {user?.email || 'guest@example.com'}
            </Typography>
          </Box>
          <Divider />
          <MenuItem onClick={() => { handleCloseUserMenu(); navigate('/profile'); }}>
            <ListItemIcon>
              <User size={18} />
            </ListItemIcon>
            <ListItemText>Profile</ListItemText>
          </MenuItem>
          <MenuItem onClick={() => { handleCloseUserMenu(); navigate('/settings'); }}>
            <ListItemIcon>
              <Settings size={18} />
            </ListItemIcon>
            <ListItemText>Settings</ListItemText>
          </MenuItem>
          <MenuItem onClick={() => { handleCloseUserMenu(); navigate('/reports'); }}>
            <ListItemIcon>
              <FileText size={18} />
            </ListItemIcon>
            <ListItemText>My Reports</ListItemText>
          </MenuItem>
          <MenuItem onClick={() => { handleCloseUserMenu(); navigate('/downloads'); }}>
            <ListItemIcon>
              <Download size={18} />
            </ListItemIcon>
            <ListItemText>Downloads</ListItemText>
          </MenuItem>
          <Divider />
          <MenuItem onClick={() => { handleCloseUserMenu(); navigate('/help'); }}>
            <ListItemIcon>
              <HelpCircle size={18} />
            </ListItemIcon>
            <ListItemText>Help & Support</ListItemText>
          </MenuItem>
          <MenuItem onClick={() => { handleCloseUserMenu(); onLogout(); }}>
            <ListItemIcon>
              <LogOut size={18} />
            </ListItemIcon>
            <ListItemText>Logout</ListItemText>
          </MenuItem>
        </Menu>
      </Toolbar>
    </AppBar>
  );
};

export default Header;
