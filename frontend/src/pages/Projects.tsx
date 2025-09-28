// Projects.tsx
import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Projects: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Projects
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Projects management page</Typography>
      </Paper>
    </Box>
  );
};

export default Projects;
