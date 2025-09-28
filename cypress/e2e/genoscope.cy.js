/**
 * E2E Tests for GenoScope Application - Continued
 */

  describe('Performance Tests', () => {
    it('should load dashboard within 3 seconds', () => {
      cy.visit('http://localhost:3000', {
        onBeforeLoad: (win) => {
          win.performance.mark('start')
        },
        onLoad: (win) => {
          win.performance.mark('end')
          win.performance.measure('pageLoad', 'start', 'end')
          const measure = win.performance.getEntriesByName('pageLoad')[0]
          expect(measure.duration).to.be.lessThan(3000)
        }
      })
    })

    it('should handle large file uploads', () => {
      cy.visit('http://localhost:3000/genomics')
      cy.contains('VCF Analysis').click()
      
      // Create a large file (mock)
      const largeFile = new File(['x'.repeat(50 * 1024 * 1024)], 'large.vcf', {
        type: 'text/plain'
      })
      
      cy.get('input[type="file"]').attachFile({
        fileContent: largeFile,
        fileName: 'large.vcf',
        mimeType: 'text/plain'
      })
      
      // Should show progress bar
      cy.get('[data-testid="upload-progress"]').should('be.visible')
      
      // Should complete eventually
      cy.contains('Upload complete', { timeout: 60000 }).should('be.visible')
    })
  })

  describe('Accessibility Tests', () => {
    it('should be navigable with keyboard', () => {
      cy.visit('http://localhost:3000')
      
      // Tab through main navigation
      cy.get('body').tab()
      cy.focused().should('have.attr', 'data-testid', 'nav-dashboard')
      
      cy.get('body').tab()
      cy.focused().should('have.attr', 'data-testid', 'nav-genomics')
      
      // Enter key should activate link
      cy.focused().type('{enter}')
      cy.url().should('include', '/genomics')
    })

    it('should have proper ARIA labels', () => {
      cy.visit('http://localhost:3000')
      
      // Check main navigation
      cy.get('nav').should('have.attr', 'aria-label', 'Main navigation')
      
      // Check buttons
      cy.get('button').each(($button) => {
        cy.wrap($button).should('have.attr', 'aria-label')
      })
      
      // Check form inputs
      cy.visit('http://localhost:3000/genomics')
      cy.get('input').each(($input) => {
        const id = $input.attr('id')
        cy.get(`label[for="${id}"]`).should('exist')
      })
    })

    it('should work with screen readers', () => {
      cy.visit('http://localhost:3000')
      
      // Check for skip navigation link
      cy.get('a[href="#main-content"]').should('exist')
      
      // Check heading hierarchy
      cy.get('h1').should('have.length', 1)
      cy.get('h2').should('have.length.at.least', 1)
      
      // Check alt text for images
      cy.get('img').each(($img) => {
        cy.wrap($img).should('have.attr', 'alt')
      })
    })
  })
})

describe('Security Tests', () => {
  it('should sanitize user input', () => {
    cy.visit('http://localhost:3000/genomics')
    cy.contains('ClinVar Search').click()
    
    // Try XSS attack
    const xssPayload = '<script>alert("XSS")</script>'
    cy.get('input[placeholder*="rsID"]').type(xssPayload)
    cy.contains('Search ClinVar').click()
    
    // Should not execute script
    cy.on('window:alert', () => {
      throw new Error('XSS vulnerability detected!')
    })
    
    // Should show sanitized error
    cy.contains('Invalid').should('be.visible')
  })

  it('should handle SQL injection attempts', () => {
    cy.visit('http://localhost:3000/genomics')
    
    const sqlPayload = "'; DROP TABLE users; --"
    cy.get('input[placeholder*="rsID"]').type(sqlPayload)
    cy.contains('Search ClinVar').click()
    
    // Should handle gracefully
    cy.contains('Invalid format').should('be.visible')
  })

  it('should enforce authentication for protected routes', () => {
    // Try to access admin panel
    cy.visit('http://localhost:3000/admin', { failOnStatusCode: false })
    
    // Should redirect to login
    cy.url().should('include', '/login')
  })

  it('should expire sessions after inactivity', () => {
    // Login first
    cy.visit('http://localhost:3000/login')
    cy.get('input[name="email"]').type('test@example.com')
    cy.get('input[name="password"]').type('password123')
    cy.contains('Login').click()
    
    // Wait for session timeout (mock with intercept)
    cy.intercept('GET', '/api/user/profile', { statusCode: 401 })
    
    // Try to perform action after timeout
    cy.visit('http://localhost:3000/genomics')
    cy.contains('VCF Analysis').click()
    
    // Should redirect to login
    cy.url().should('include', '/login')
    cy.contains('Session expired').should('be.visible')
  })
})
