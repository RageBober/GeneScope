#!/bin/bash

echo "🚀 Starting GenoScope monitoring stack..."

# Start monitoring services
cd monitoring
docker-compose up -d

echo "✅ Monitoring services started:"
echo "📊 Grafana:     http://localhost:3001 (admin/admin)"
echo "📈 Prometheus:  http://localhost:9090"
echo "📋 AlertManager: http://localhost:9093"
echo "📊 Node Metrics: http://localhost:9100"
echo "🐳 cAdvisor:     http://localhost:8080"

echo ""
echo "💡 To integrate with main application:"
echo "   Add 'monitoring' network to your main docker-compose.yml"
echo "   or use docker-compose.full.yml"
