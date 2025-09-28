import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Search: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Search
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Search page</Typography>
      </Paper>
    </Box>
  );
};

export default Search;
