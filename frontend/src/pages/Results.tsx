import React from 'react';
import { Box, Typography, Paper } from '@mui/material';
import { useParams } from 'react-router-dom';

const Results: React.FC = () => {
  const { id } = useParams();

  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Analysis Results #{id}
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Results page content will be here</Typography>
      </Paper>
    </Box>
  );
};

export default Results;
