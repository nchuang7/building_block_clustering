
# GCS bucket and file paths
GCS_BUCKET = "satomic-building-blocks"
GCS_BASE_PATH = f"gs://{GCS_BUCKET}"

ACIDS_PATH = f'{GCS_BASE_PATH}/acids.parquet'
AMINES_PATH = f'{GCS_BASE_PATH}/amines.parquet'

# app settings
#CACHE_TTL = 86400 # cache for 24 hours
