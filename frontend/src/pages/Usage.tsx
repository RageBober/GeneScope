import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Usage: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Usage Analytics
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Usage analytics page</Typography>
      </Paper>
    </Box>
  );
};

export default Usage;
