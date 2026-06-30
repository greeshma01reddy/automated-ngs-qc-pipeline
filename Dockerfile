# Use a lightweight Python base image
FROM python:3.10-slim

# Install system dependencies including Java 21 (for Nextflow) and curl
RUN apt-get update && apt-get install -y \
    openjdk-21-jre-headless \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow globally
RUN curl -s https://get.nextflow.io | bash && chmod +x nextflow && mv nextflow /usr/local/bin/

# Set the working directory inside the container
WORKDIR /app

# Copy dependency mappings first to optimize build caching
COPY requirements.txt .

# Install Python scientific stack & Streamlit
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your pipeline architecture code
COPY . .

# Expose Streamlit's default production port
EXPOSE 8501

# Command to spin up the dashboard by default
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
