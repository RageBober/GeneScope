# GenoScope - Genomics Analysis Platform

## 🚀 Quick Start

### Prerequisites
- Python 3.8+
- Docker & Docker Compose
- Node.js 16+ (for frontend development)

### Development Setup

1. **Clone and setup environment:**
   ```bash
   git clone <repository>
   cd BioForge_edit_branch
   python -m venv .venv
   source .venv/bin/activate  # Linux/Mac
   # or
   .venv\Scripts\activate  # Windows
   ```

2. **Install dependencies:**
   ```bash
   pip install -e .
   # or
   pip install -r requirements.txt
   ```

3. **Start services (Docker):**
   ```bash
   # For full stack
   docker-compose up -d
   
   # For development (minimal)
   docker-compose -f docker-compose.dev.yml up -d
   ```

4. **Start frontend:**
   ```bash
   cd frontend
   npm install
   npm start
   ```

### Project Structure

```
BioForge_edit_branch/
├── src/                    # Python source code
├── frontend/              # React frontend
├── data/
│   ├── references/        # Reference genomes
│   ├── samples/          # Sample data
│   └── reports/          # Analysis reports
├── tests/                # Test files
├── docs/                 # Documentation
├── scripts/              # Utility scripts
└── config/               # Configuration files
```

### Configuration

- **Development:** Use `docker-compose.dev.yml` for minimal setup
- **Production:** Use `docker-compose.yml` for full stack
- **Local:** Modify `.env` to use SQLite for local development

### Database Migration

```bash
# If using PostgreSQL
docker-compose exec backend alembic upgrade head

# If using SQLite (local)
python -m alembic upgrade head
```

### Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src --cov-report=html

# Run specific test types
pytest -m "not slow"  # Skip slow tests
pytest -m integration  # Only integration tests
```

### Reference Data

Large reference files (like human genome) should be downloaded separately:

```bash
# Download to correct location
wget -P data/references/human/ https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### Useful Commands

- `./quick_start.sh` - Start development environment
- `./start.sh` - Start backend server
- `./stop.sh` - Stop all services
- `./run_tests.py` - Run test suite

## 🔧 Development

### Code Style
- Black for formatting
- Flake8 for linting  
- MyPy for type checking
- Pytest for testing

### Pre-commit Hooks
```bash
pip install pre-commit
pre-commit install
```

## 📊 Monitoring

- **API Docs:** http://localhost:8000/docs
- **Frontend:** http://localhost:3000
- **Flower (Celery):** http://localhost:5555 (full stack only)
- **Database:** localhost:5432

## 🤝 Contributing

1. Create feature branch
2. Make changes
3. Run tests: `pytest`
4. Submit pull request

## 📝 License

See LICENSE file for details.
