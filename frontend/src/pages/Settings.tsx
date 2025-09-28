import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Settings: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Settings
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Settings page content will be here</Typography>
      </Paper>
    </Box>
  );
};

export default Settings;
