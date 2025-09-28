import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const QuickRun: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Quick Run
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Quick Run analysis page</Typography>
      </Paper>
    </Box>
  );
};

export default QuickRun;
