import React, { useState, useEffect, createContext } from 'react';
import { Routes, Route, Navigate, useLocation } from 'react-router-dom';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import { CssBaseline, Box } from '@mui/material';
import { SnackbarProvider } from 'notistack';

// Import components
import Header from './components/Layout/Header';
import Sidebar from './components/Layout/Sidebar';
import Login from './pages/Login';
import Dashboard from './pages/Dashboard';
import Analysis from './pages/Analysis';
import Results from './pages/Results';
import Reports from './pages/Reports';
import Projects from './pages/Projects';
import Samples from './pages/Samples';
import Settings from './pages/Settings';
import Billing from './pages/Billing';
import QuickRun from './pages/QuickRun';
import Pipeline from './pages/Pipeline';
import Compare from './pages/Compare';
import Search from './pages/Search';
import Storage from './pages/Storage';
import Usage from './pages/Usage';
import Datasets from './pages/Datasets';
import References from './pages/References';
import GenomicsAnalysis from './pages/GenomicsAnalysis';

// Create theme
const theme = createTheme({
  palette: {
    mode: 'light',
    primary: {
      main: '#4F46E5',
      light: '#818CF8',
      dark: '#3730A3',
    },
    secondary: {
      main: '#EC4899',
      light: '#F9A8D4',
      dark: '#BE185D',
    },
    success: {
      main: '#10B981',
    },
    error: {
      main: '#EF4444',
    },
    warning: {
      main: '#F59E0B',
    },
    info: {
      main: '#3B82F6',
    },
    background: {
      default: '#F9FAFB',
      paper: '#FFFFFF',
    },
  },
  typography: {
    fontFamily: '"Inter", "Roboto", "Helvetica", "Arial", sans-serif',
    h1: {
      fontSize: '2.5rem',
      fontWeight: 700,
    },
    h2: {
      fontSize: '2rem',
      fontWeight: 600,
    },
    h3: {
      fontSize: '1.75rem',
      fontWeight: 600,
    },
    h4: {
      fontSize: '1.5rem',
      fontWeight: 600,
    },
    h5: {
      fontSize: '1.25rem',
      fontWeight: 600,
    },
    button: {
      textTransform: 'none',
    },
  },
  shape: {
    borderRadius: 8,
  },
});

// User Context
export const UserContext = createContext<any>(null);

function App() {
  const [isAuthenticated, setIsAuthenticated] = useState(false);
  const [user, setUser] = useState<any>(null);
  const [sidebarOpen, setSidebarOpen] = useState(true);
  const [loading, setLoading] = useState(true);
  const location = useLocation();

  useEffect(() => {
    // Check authentication status
    const token = localStorage.getItem('access_token');
    if (token) {
      setIsAuthenticated(true);
      // Mock user data - replace with actual API call
      setUser({
        name: 'John Doe',
        email: 'john.doe@example.com',
        role: 'researcher',
      });
    }
    setLoading(false);
  }, []);

  const handleLogout = () => {
    localStorage.removeItem('access_token');
    localStorage.removeItem('refresh_token');
    setIsAuthenticated(false);
    setUser(null);
  };

  const handleLogin = (token: string, userData: any) => {
    localStorage.setItem('access_token', token);
    setIsAuthenticated(true);
    setUser(userData);
  };

  if (loading) {
    return (
      <ThemeProvider theme={theme}>
        <CssBaseline />
        <Box
          sx={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            height: '100vh',
          }}
        >
          Loading...
        </Box>
      </ThemeProvider>
    );
  }

  const isLoginPage = location.pathname === '/login';

  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <SnackbarProvider maxSnack={3}>
        <UserContext.Provider value={{ user, setUser }}>
          <Box sx={{ display: 'flex', minHeight: '100vh' }}>
            {isAuthenticated && !isLoginPage && (
              <>
                <Header
                  onMenuClick={() => setSidebarOpen(!sidebarOpen)}
                  onLogout={handleLogout}
                />
                <Sidebar open={sidebarOpen} />
              </>
            )}
            <Box
              component="main"
              sx={{
                flexGrow: 1,
                pt: isAuthenticated && !isLoginPage ? 8 : 0,
                pl: isAuthenticated && !isLoginPage && sidebarOpen ? '240px' : 0,
                transition: 'padding-left 0.3s',
                backgroundColor: theme.palette.background.default,
                minHeight: '100vh',
              }}
            >
              <Box sx={{ p: isAuthenticated && !isLoginPage ? 3 : 0 }}>
                <Routes>
                  <Route
                    path="/login"
                    element={
                      isAuthenticated ? (
                        <Navigate to="/dashboard" />
                      ) : (
                        <Login onLogin={handleLogin} />
                      )
                    }
                  />
                  <Route
                    path="/dashboard"
                    element={
                      isAuthenticated ? <Dashboard /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/analysis"
                    element={
                      isAuthenticated ? <Analysis /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/results/:id"
                    element={
                      isAuthenticated ? <Results /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/reports"
                    element={
                      isAuthenticated ? <Reports /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/projects"
                    element={
                      isAuthenticated ? <Projects /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/samples"
                    element={
                      isAuthenticated ? <Samples /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/quick-run"
                    element={
                      isAuthenticated ? <QuickRun /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/pipeline"
                    element={
                      isAuthenticated ? <Pipeline /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/compare"
                    element={
                      isAuthenticated ? <Compare /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/search"
                    element={
                      isAuthenticated ? <Search /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/storage"
                    element={
                      isAuthenticated ? <Storage /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/usage"
                    element={
                      isAuthenticated ? <Usage /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/datasets"
                    element={
                      isAuthenticated ? <Datasets /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/references"
                    element={
                      isAuthenticated ? <References /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/settings"
                    element={
                      isAuthenticated ? <Settings /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/billing"
                    element={
                      isAuthenticated ? <Billing /> : <Navigate to="/login" />
                    }
                  />
                  <Route
                    path="/genomics"
                    element={
                      isAuthenticated ? <GenomicsAnalysis /> : <Navigate to="/login" />
                    }
                  />
                  <Route path="/" element={<Navigate to="/dashboard" />} />
                </Routes>
              </Box>
            </Box>
          </Box>
        </UserContext.Provider>
      </SnackbarProvider>
    </ThemeProvider>
  );
}

export default App;
