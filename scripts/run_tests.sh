#!/bin/bash

# GenoScope Comprehensive Test Suite Runner
# Запускает все типы тестов и генерирует отчеты

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test results
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_header() {
    echo ""
    echo "=================================="
    echo "$1"
    echo "=================================="
    echo ""
}

# Parse arguments
RUN_UNIT=true
RUN_INTEGRATION=true
RUN_E2E=true
RUN_PERFORMANCE=false
RUN_SECURITY=false
COVERAGE_REPORT=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --unit-only)
            RUN_INTEGRATION=false
            RUN_E2E=false
            shift
            ;;
        --integration-only)
            RUN_UNIT=false
            RUN_E2E=false
            shift
            ;;
        --e2e-only)
            RUN_UNIT=false
            RUN_INTEGRATION=false
            shift
            ;;
        --performance)
            RUN_PERFORMANCE=true
            shift
            ;;
        --security)
            RUN_SECURITY=true
            shift
            ;;
        --all)
            RUN_PERFORMANCE=true
            RUN_SECURITY=true
            shift
            ;;
        --no-coverage)
            COVERAGE_REPORT=false
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--unit-only|--integration-only|--e2e-only|--performance|--security|--all|--no-coverage]"
            exit 1
            ;;
    esac
done

# Start time
START_TIME=$(date +%s)

# Create test results directory
mkdir -p test-results

print_header "GenoScope Test Suite"
log_info "Starting comprehensive test suite..."

# 1. Unit Tests
if [ "$RUN_UNIT" = true ]; then
    print_header "Running Unit Tests"
    
    if pytest tests/unit/ -v --tb=short --junit-xml=test-results/unit-tests.xml; then
        log_success "Unit tests passed!"
        ((PASSED_TESTS++))
    else
        log_error "Unit tests failed!"
        ((FAILED_TESTS++))
    fi
    ((TOTAL_TESTS++))
fi

# 2. Integration Tests
if [ "$RUN_INTEGRATION" = true ]; then
    print_header "Running Integration Tests"
    
    # Start test database
    log_info "Starting test database..."
    docker-compose -f docker-compose.test.yml up -d postgres redis
    sleep 5
    
    if pytest tests/integration/ -v --tb=short --junit-xml=test-results/integration-tests.xml; then
        log_success "Integration tests passed!"
        ((PASSED_TESTS++))
    else
        log_error "Integration tests failed!"
        ((FAILED_TESTS++))
    fi
    ((TOTAL_TESTS++))
    
    # Stop test database
    docker-compose -f docker-compose.test.yml down
fi

# 3. E2E Tests
if [ "$RUN_E2E" = true ]; then
    print_header "Running E2E Tests"
    
    # Start application
    log_info "Starting application for E2E tests..."
    docker-compose up -d
    
    # Wait for services to be ready
    log_info "Waiting for services to be ready..."
    for i in {1..30}; do
        if curl -f http://localhost:8000/health > /dev/null 2>&1; then
            break
        fi
        sleep 2
    done
    
    # Run Cypress tests
    cd cypress
    npm install
    
    if npm run cy:run; then
        log_success "E2E tests passed!"
        ((PASSED_TESTS++))
    else
        log_error "E2E tests failed!"
        ((FAILED_TESTS++))
    fi
    ((TOTAL_TESTS++))
    
    cd ..
    docker-compose down
fi

# 4. Performance Tests
if [ "$RUN_PERFORMANCE" = true ]; then
    print_header "Running Performance Tests"
    
    # Start application for performance testing
    log_info "Starting application for performance tests..."
    docker-compose up -d
    sleep 10
    
    # Run Locust tests
    log_info "Running load tests with Locust..."
    
    if locust -f tests/performance/locustfile.py \
              --host=http://localhost:8000 \
              --users=50 \
              --spawn-rate=5 \
              --run-time=60s \
              --headless \
              --html=test-results/performance-report.html \
              --csv=test-results/performance; then
        
        # Check performance criteria
        AVG_RESPONSE=$(grep "Aggregated" test-results/performance_stats.csv | cut -d',' -f6)
        if (( $(echo "$AVG_RESPONSE < 500" | bc -l) )); then
            log_success "Performance tests passed! Avg response: ${AVG_RESPONSE}ms"
            ((PASSED_TESTS++))
        else
            log_error "Performance tests failed! Avg response: ${AVG_RESPONSE}ms (SLA: <500ms)"
            ((FAILED_TESTS++))
        fi
    else
        log_error "Performance tests failed to run!"
        ((FAILED_TESTS++))
    fi
    ((TOTAL_TESTS++))
    
    docker-compose down
fi

# 5. Security Tests
if [ "$RUN_SECURITY" = true ]; then
    print_header "Running Security Tests"
    
    # Run Bandit for Python security
    log_info "Running Bandit security scan..."
    if bandit -r src/ -f json -o test-results/bandit-report.json; then
        log_success "No security issues found by Bandit"
    else
        log_warning "Bandit found potential security issues"
    fi
    
    # Run Safety for dependency check
    log_info "Checking dependencies for vulnerabilities..."
    if safety check --json > test-results/safety-report.json; then
        log_success "No vulnerable dependencies found"
    else
        log_warning "Vulnerable dependencies detected"
    fi
    
    # Run OWASP ZAP if available
    if command -v zap-cli &> /dev/null; then
        log_info "Running OWASP ZAP security scan..."
        docker-compose up -d
        sleep 10
        
        zap-cli start --start-options '-config api.disablekey=true'
        zap-cli open-url http://localhost:8000
        zap-cli spider http://localhost:8000
        zap-cli active-scan http://localhost:8000
        zap-cli report -o test-results/zap-report.html -f html
        
        ALERTS=$(zap-cli alerts -l High)
        if [ -z "$ALERTS" ]; then
            log_success "No high-risk vulnerabilities found"
            ((PASSED_TESTS++))
        else
            log_error "High-risk vulnerabilities detected!"
            ((FAILED_TESTS++))
        fi
        
        zap-cli shutdown
        docker-compose down
    else
        log_warning "OWASP ZAP not installed, skipping web security scan"
    fi
    ((TOTAL_TESTS++))
fi

# 6. Coverage Report
if [ "$COVERAGE_REPORT" = true ]; then
    print_header "Generating Coverage Report"
    
    # Combine coverage from all test runs
    coverage combine
    coverage report
    coverage html
    
    COVERAGE=$(coverage report | grep TOTAL | awk '{print $4}' | sed 's/%//')
    
    if (( $(echo "$COVERAGE >= 80" | bc -l) )); then
        log_success "Code coverage: ${COVERAGE}% (Target: ≥80%)"
    else
        log_warning "Code coverage: ${COVERAGE}% (Target: ≥80%)"
    fi
    
    log_info "HTML coverage report generated in htmlcov/"
fi

# 7. Generate consolidated report
print_header "Test Results Summary"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

cat > test-results/summary.json <<EOF
{
  "timestamp": "$(date -Iseconds)",
  "duration_seconds": $DURATION,
  "total_suites": $TOTAL_TESTS,
  "passed_suites": $PASSED_TESTS,
  "failed_suites": $FAILED_TESTS,
  "coverage_percentage": ${COVERAGE:-0},
  "unit_tests": "$([[ $RUN_UNIT == true ]] && echo "executed" || echo "skipped")",
  "integration_tests": "$([[ $RUN_INTEGRATION == true ]] && echo "executed" || echo "skipped")",
  "e2e_tests": "$([[ $RUN_E2E == true ]] && echo "executed" || echo "skipped")",
  "performance_tests": "$([[ $RUN_PERFORMANCE == true ]] && echo "executed" || echo "skipped")",
  "security_tests": "$([[ $RUN_SECURITY == true ]] && echo "executed" || echo "skipped")"
}
EOF

# Print summary
echo ""
echo "========================================="
echo "            TEST SUMMARY                 "
echo "========================================="
echo ""
echo "Duration: ${DURATION} seconds"
echo "Total Test Suites: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed: ${RED}$FAILED_TESTS${NC}"
if [ "$COVERAGE_REPORT" = true ]; then
    echo "Code Coverage: ${COVERAGE}%"
fi
echo ""

# Generate HTML report
cat > test-results/index.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>GenoScope Test Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }
        .summary { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }
        .metric { background: #f8f9fa; padding: 20px; border-radius: 5px; text-align: center; }
        .metric.passed { border-left: 5px solid #28a745; }
        .metric.failed { border-left: 5px solid #dc3545; }
        .metric h3 { margin: 0; color: #495057; }
        .metric .value { font-size: 2em; font-weight: bold; margin: 10px 0; }
        .passed .value { color: #28a745; }
        .failed .value { color: #dc3545; }
        .reports { margin: 30px 0; }
        .reports a { display: inline-block; margin: 10px; padding: 10px 20px; background: #007bff; color: white; text-decoration: none; border-radius: 5px; }
    </style>
</head>
<body>
    <div class="header">
        <h1>GenoScope Test Report</h1>
        <p>Generated: $(date)</p>
    </div>
    
    <div class="summary">
        <div class="metric passed">
            <h3>Passed</h3>
            <div class="value">$PASSED_TESTS</div>
        </div>
        <div class="metric failed">
            <h3>Failed</h3>
            <div class="value">$FAILED_TESTS</div>
        </div>
        <div class="metric">
            <h3>Coverage</h3>
            <div class="value">${COVERAGE:-N/A}%</div>
        </div>
        <div class="metric">
            <h3>Duration</h3>
            <div class="value">${DURATION}s</div>
        </div>
    </div>
    
    <div class="reports">
        <h2>Detailed Reports</h2>
        <a href="../htmlcov/index.html">Coverage Report</a>
        <a href="unit-tests.xml">Unit Test Results</a>
        <a href="integration-tests.xml">Integration Test Results</a>
        <a href="performance-report.html">Performance Report</a>
        <a href="summary.json">JSON Summary</a>
    </div>
</body>
</html>
EOF

log_success "Test report generated: test-results/index.html"

# Exit with appropriate code
if [ $FAILED_TESTS -gt 0 ]; then
    log_error "Test suite failed with $FAILED_TESTS failures"
    exit 1
else
    log_success "All tests passed successfully!"
    exit 0
fi
