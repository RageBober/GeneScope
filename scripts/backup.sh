#!/bin/bash

# GenoScope Automated Backup Script
# Runs database backups, file system backups, and uploads to S3

set -e

# Configuration
BACKUP_DIR="./data/backups"
S3_BUCKET="s3://genoscope-backups"
DATE=$(date +%Y%m%d_%H%M%S)
RETENTION_DAYS=30
SLACK_WEBHOOK_URL="${SLACK_WEBHOOK_URL}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1"
}

warning() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1"
}

# Send notification to Slack
notify_slack() {
    local message=$1
    local status=$2
    
    if [ ! -z "$SLACK_WEBHOOK_URL" ]; then
        curl -X POST -H 'Content-type: application/json' \
            --data "{\"text\":\"GenoScope Backup: $status - $message\"}" \
            $SLACK_WEBHOOK_URL
    fi
}

# Create backup directory
mkdir -p $BACKUP_DIR
cd $BACKUP_DIR

log "Starting GenoScope backup process..."

# 1. Database Backup
backup_database() {
    log "Backing up PostgreSQL database..."
    
    DB_BACKUP_FILE="genoscope_db_${DATE}.sql.gz"
    
    # Dump database
    PGPASSWORD=$DB_PASSWORD pg_dump \
        -h $DB_HOST \
        -U $DB_USER \
        -d $DB_NAME \
        --no-owner \
        --no-acl \
        --clean \
        --if-exists | gzip > $DB_BACKUP_FILE
    
    if [ $? -eq 0 ]; then
        log "Database backup completed: $DB_BACKUP_FILE"
        
        # Upload to S3
        aws s3 cp $DB_BACKUP_FILE $S3_BUCKET/database/ \
            --storage-class STANDARD_IA \
            --server-side-encryption AES256
        
        if [ $? -eq 0 ]; then
            log "Database backup uploaded to S3"
            rm $DB_BACKUP_FILE
        else
            error "Failed to upload database backup to S3"
            return 1
        fi
    else
        error "Database backup failed"
        return 1
    fi
}

# 2. Application Data Backup
backup_application_data() {
    log "Backing up application data..."
    
    DATA_BACKUP_FILE="genoscope_data_${DATE}.tar.gz"
    
    # Create tar archive of important directories
    tar -czf $DATA_BACKUP_FILE \
        /data/pipeline_runs \
        /data/uploads \
        /data/reports \
        --exclude='*.tmp' \
        --exclude='*/cache/*'
    
    if [ $? -eq 0 ]; then
        log "Application data backup completed: $DATA_BACKUP_FILE"
        
        # Upload to S3
        aws s3 cp $DATA_BACKUP_FILE $S3_BUCKET/application-data/ \
            --storage-class STANDARD_IA \
            --server-side-encryption AES256
        
        if [ $? -eq 0 ]; then
            log "Application data backup uploaded to S3"
            rm $DATA_BACKUP_FILE
        else
            error "Failed to upload application data backup to S3"
            return 1
        fi
    else
        error "Application data backup failed"
        return 1
    fi
}

# 3. Configuration Backup
backup_configuration() {
    log "Backing up configuration files..."
    
    CONFIG_BACKUP_FILE="genoscope_config_${DATE}.tar.gz"
    
    # Backup configuration files
    tar -czf $CONFIG_BACKUP_FILE \
        /app/.env \
        /app/config \
        /etc/nginx/sites-available/genoscope \
        /etc/systemd/system/genoscope*.service \
        --exclude='*secret*' \
        --exclude='*key*'
    
    if [ $? -eq 0 ]; then
        log "Configuration backup completed: $CONFIG_BACKUP_FILE"
        
        # Encrypt sensitive configuration
        gpg --encrypt --recipient backup@genoscope.com $CONFIG_BACKUP_FILE
        
        # Upload to S3
        aws s3 cp ${CONFIG_BACKUP_FILE}.gpg $S3_BUCKET/configuration/ \
            --storage-class STANDARD_IA \
            --server-side-encryption AES256
        
        if [ $? -eq 0 ]; then
            log "Configuration backup uploaded to S3"
            rm $CONFIG_BACKUP_FILE ${CONFIG_BACKUP_FILE}.gpg
        else
            error "Failed to upload configuration backup to S3"
            return 1
        fi
    else
        error "Configuration backup failed"
        return 1
    fi
}

# 4. Redis Backup
backup_redis() {
    log "Backing up Redis data..."
    
    REDIS_BACKUP_FILE="genoscope_redis_${DATE}.rdb"
    
    # Trigger Redis BGSAVE
    redis-cli -h $REDIS_HOST BGSAVE
    
    # Wait for backup to complete
    while [ $(redis-cli -h $REDIS_HOST LASTSAVE) -eq $(redis-cli -h $REDIS_HOST LASTSAVE) ]; do
        sleep 1
    done
    
    # Copy dump file
    cp /var/lib/redis/dump.rdb $REDIS_BACKUP_FILE
    
    if [ $? -eq 0 ]; then
        log "Redis backup completed: $REDIS_BACKUP_FILE"
        
        # Compress and upload to S3
        gzip $REDIS_BACKUP_FILE
        aws s3 cp ${REDIS_BACKUP_FILE}.gz $S3_BUCKET/redis/ \
            --storage-class STANDARD_IA \
            --server-side-encryption AES256
        
        if [ $? -eq 0 ]; then
            log "Redis backup uploaded to S3"
            rm ${REDIS_BACKUP_FILE}.gz
        else
            error "Failed to upload Redis backup to S3"
            return 1
        fi
    else
        error "Redis backup failed"
        return 1
    fi
}

# 5. Cleanup old backups
cleanup_old_backups() {
    log "Cleaning up old backups (retention: $RETENTION_DAYS days)..."
    
    # Calculate cutoff date
    CUTOFF_DATE=$(date -d "$RETENTION_DAYS days ago" +%s)
    
    # List and delete old backups from S3
    aws s3 ls $S3_BUCKET/ --recursive | while read -r line; do
        CREATE_DATE=$(echo $line | awk '{print $1" "$2}')
        CREATE_DATE_EPOCH=$(date -d "$CREATE_DATE" +%s)
        FILE_PATH=$(echo $line | awk '{print $4}')
        
        if [ $CREATE_DATE_EPOCH -lt $CUTOFF_DATE ]; then
            log "Deleting old backup: $FILE_PATH"
            aws s3 rm $S3_BUCKET/$FILE_PATH
        fi
    done
}

# 6. Verify backups
verify_backups() {
    log "Verifying backups..."
    
    # Check if today's backups exist in S3
    TODAY=$(date +%Y%m%d)
    
    DB_EXISTS=$(aws s3 ls $S3_BUCKET/database/ | grep $TODAY | wc -l)
    DATA_EXISTS=$(aws s3 ls $S3_BUCKET/application-data/ | grep $TODAY | wc -l)
    CONFIG_EXISTS=$(aws s3 ls $S3_BUCKET/configuration/ | grep $TODAY | wc -l)
    REDIS_EXISTS=$(aws s3 ls $S3_BUCKET/redis/ | grep $TODAY | wc -l)
    
    if [ $DB_EXISTS -gt 0 ] && [ $DATA_EXISTS -gt 0 ] && [ $CONFIG_EXISTS -gt 0 ] && [ $REDIS_EXISTS -gt 0 ]; then
        log "All backups verified successfully"
        return 0
    else
        error "Backup verification failed"
        return 1
    fi
}

# 7. Generate backup report
generate_report() {
    log "Generating backup report..."
    
    REPORT_FILE="backup_report_${DATE}.txt"
    
    cat > $REPORT_FILE << EOF
GenoScope Backup Report
========================
Date: $(date)
Status: $1

Backup Summary:
--------------
Database Backup: $2
Application Data: $3
Configuration: $4
Redis: $5

S3 Bucket Usage:
---------------
$(aws s3 ls $S3_BUCKET/ --recursive --summarize | tail -2)

Recent Backups:
--------------
$(aws s3 ls $S3_BUCKET/ --recursive | tail -20)

EOF
    
    # Email report
    if [ ! -z "$ADMIN_EMAIL" ]; then
        mail -s "GenoScope Backup Report - $(date +%Y-%m-%d)" $ADMIN_EMAIL < $REPORT_FILE
    fi
    
    # Upload report to S3
    aws s3 cp $REPORT_FILE $S3_BUCKET/reports/
    
    rm $REPORT_FILE
}

# Main execution
main() {
    local overall_status="SUCCESS"
    local db_status="✓"
    local data_status="✓"
    local config_status="✓"
    local redis_status="✓"
    
    # Run backups
    if ! backup_database; then
        db_status="✗"
        overall_status="PARTIAL"
    fi
    
    if ! backup_application_data; then
        data_status="✗"
        overall_status="PARTIAL"
    fi
    
    if ! backup_configuration; then
        config_status="✗"
        overall_status="PARTIAL"
    fi
    
    if ! backup_redis; then
        redis_status="✗"
        overall_status="PARTIAL"
    fi
    
    # Cleanup and verify
    cleanup_old_backups
    
    if ! verify_backups; then
        overall_status="FAILED"
    fi
    
    # Generate report
    generate_report "$overall_status" "$db_status" "$data_status" "$config_status" "$redis_status"
    
    # Send notification
    notify_slack "Backup completed with status: $overall_status" "$overall_status"
    
    if [ "$overall_status" == "SUCCESS" ]; then
        log "Backup process completed successfully"
        exit 0
    else
        error "Backup process completed with errors"
        exit 1
    fi
}

# Run main function
main
