import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Storage: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Storage
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Storage management page</Typography>
      </Paper>
    </Box>
  );
};

export default Storage;
