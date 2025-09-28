import React from 'react';
import {
  Drawer,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  ListItemButton,
  Divider,
  Typography,
  Box,
  Collapse,
  Chip,
} from '@mui/material';
import { useNavigate, useLocation } from 'react-router-dom';
import {
  Home,
  Activity,
  FileText,
  FolderOpen,
  TestTube2,
  Settings,
  CreditCard,
  PlayCircle,
  GitBranch,
  Search,
  BarChart3,
  HardDrive,
  TrendingUp,
  Database,
  Dna,
  ChevronDown,
  ChevronRight,
  Beaker,
  Package,
} from 'lucide-react';

interface SidebarProps {
  open: boolean;
}

const Sidebar: React.FC<SidebarProps> = ({ open }) => {
  const navigate = useNavigate();
  const location = useLocation();
  const [expandedItems, setExpandedItems] = React.useState<string[]>(['analysis']);

  const handleToggle = (item: string) => {
    setExpandedItems((prev) =>
      prev.includes(item)
        ? prev.filter((i) => i !== item)
        : [...prev, item]
    );
  };

  const menuItems = [
    {
      id: 'dashboard',
      title: 'Dashboard',
      icon: <Home size={20} />,
      path: '/dashboard',
    },
    {
      id: 'analysis',
      title: 'Analysis',
      icon: <Activity size={20} />,
      children: [
        {
          id: 'new-analysis',
          title: 'New Analysis',
          path: '/analysis',
          badge: 'New',
        },
        {
          id: 'quick-run',
          title: 'Quick Run',
          path: '/quick-run',
          icon: <PlayCircle size={18} />,
        },
        {
          id: 'pipeline',
          title: 'Pipeline Builder',
          path: '/pipeline',
          icon: <GitBranch size={18} />,
        },
        {
          id: 'compare',
          title: 'Compare Analyses',
          path: '/compare',
          icon: <BarChart3 size={18} />,
        },
      ],
    },
    {
      id: 'results',
      title: 'Results',
      icon: <FileText size={20} />,
      path: '/reports',
      badge: '3',
    },
    {
      id: 'data',
      title: 'Data Management',
      icon: <Database size={20} />,
      children: [
        {
          id: 'projects',
          title: 'Projects',
          path: '/projects',
          icon: <FolderOpen size={18} />,
        },
        {
          id: 'samples',
          title: 'Samples',
          path: '/samples',
          icon: <TestTube2 size={18} />,
        },
        {
          id: 'datasets',
          title: 'Datasets',
          path: '/datasets',
          icon: <Package size={18} />,
        },
        {
          id: 'references',
          title: 'Reference Genomes',
          path: '/references',
          icon: <Dna size={18} />,
        },
      ],
    },
    {
      id: 'genomics',
      title: 'Genomics Analysis',
      icon: <Beaker size={20} />,
      path: '/genomics',
      badge: 'New',
    },
    {
      id: 'search',
      title: 'Search',
      icon: <Search size={20} />,
      path: '/search',
    },
    {
      id: 'resources',
      title: 'Resources',
      icon: <TrendingUp size={20} />,
      children: [
        {
          id: 'storage',
          title: 'Storage',
          path: '/storage',
          icon: <HardDrive size={18} />,
        },
        {
          id: 'usage',
          title: 'Usage Analytics',
          path: '/usage',
          icon: <BarChart3 size={18} />,
        },
      ],
    },
  ];

  const bottomMenuItems = [
    {
      id: 'settings',
      title: 'Settings',
      icon: <Settings size={20} />,
      path: '/settings',
    },
    {
      id: 'billing',
      title: 'Billing',
      icon: <CreditCard size={20} />,
      path: '/billing',
    },
  ];

  const drawerWidth = 240;

  const renderMenuItem = (item: any, nested = false) => {
    const isActive = location.pathname === item.path;
    const hasChildren = item.children && item.children.length > 0;
    const isExpanded = expandedItems.includes(item.id);

    if (hasChildren) {
      return (
        <React.Fragment key={item.id}>
          <ListItemButton
            onClick={() => handleToggle(item.id)}
            sx={{
              pl: nested ? 4 : 2,
              '&:hover': {
                backgroundColor: 'action.hover',
              },
            }}
          >
            <ListItemIcon sx={{ minWidth: 40 }}>
              {item.icon}
            </ListItemIcon>
            <ListItemText primary={item.title} />
            {isExpanded ? <ChevronDown size={18} /> : <ChevronRight size={18} />}
          </ListItemButton>
          <Collapse in={isExpanded} timeout="auto" unmountOnExit>
            <List component="div" disablePadding>
              {item.children.map((child: any) => renderMenuItem(child, true))}
            </List>
          </Collapse>
        </React.Fragment>
      );
    }

    return (
      <ListItemButton
        key={item.id}
        onClick={() => navigate(item.path)}
        sx={{
          pl: nested ? 6 : 2,
          backgroundColor: isActive ? 'action.selected' : 'transparent',
          '&:hover': {
            backgroundColor: 'action.hover',
          },
        }}
      >
        {item.icon && (
          <ListItemIcon sx={{ minWidth: 40 }}>
            {item.icon}
          </ListItemIcon>
        )}
        <ListItemText primary={item.title} />
        {item.badge && (
          <Chip
            label={item.badge}
            size="small"
            color={item.badge === 'New' ? 'primary' : 'default'}
            sx={{ height: 20, fontSize: '0.75rem' }}
          />
        )}
      </ListItemButton>
    );
  };

  return (
    <Drawer
      variant="persistent"
      anchor="left"
      open={open}
      sx={{
        width: open ? drawerWidth : 0,
        flexShrink: 0,
        '& .MuiDrawer-paper': {
          width: drawerWidth,
          boxSizing: 'border-box',
          top: '64px',
          height: 'calc(100% - 64px)',
          borderRight: '1px solid rgba(0, 0, 0, 0.12)',
        },
      }}
    >
      <Box sx={{ overflow: 'auto', height: '100%', display: 'flex', flexDirection: 'column' }}>
        {/* Main Menu Items */}
        <List sx={{ flexGrow: 1 }}>
          {menuItems.map((item) => renderMenuItem(item))}
        </List>

        {/* Bottom Menu Items */}
        <Box>
          <Divider />
          <List>
            {bottomMenuItems.map((item) => renderMenuItem(item))}
          </List>
          
          {/* Version Info */}
          <Box sx={{ p: 2, textAlign: 'center' }}>
            <Typography variant="caption" color="text.secondary">
              GenoScope v1.0.0-beta
            </Typography>
          </Box>
        </Box>
      </Box>
    </Drawer>
  );
};

export default Sidebar;
