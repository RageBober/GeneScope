const { defineConfig } = require('cypress')

module.exports = defineConfig({
  e2e: {
    baseUrl: 'http://localhost:3000',
    viewportWidth: 1280,
    viewportHeight: 720,
    video: true,
    screenshotOnRunFailure: true,
    
    setupNodeEvents(on, config) {
      // Implement node event listeners here
      
      // Code coverage
      require('@cypress/code-coverage/task')(on, config)
      
      // File upload support
      on('task', {
        // Custom tasks
        log(message) {
          console.log(message)
          return null
        },
        table(message) {
          console.table(message)
          return null
        }
      })
      
      return config
    },
    
    env: {
      apiUrl: 'http://localhost:8000',
      coverage: true
    },
    
    defaultCommandTimeout: 10000,
    requestTimeout: 10000,
    responseTimeout: 10000,
    
    retries: {
      runMode: 2,
      openMode: 0
    },
    
    experimentalStudio: true,
    experimentalWebKitSupport: true,
  },
  
  component: {
    devServer: {
      framework: 'react',
      bundler: 'webpack',
    },
    specPattern: 'src/**/*.cy.{js,jsx,ts,tsx}',
  },
})
