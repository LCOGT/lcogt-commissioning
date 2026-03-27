FROM python:3.12
# Configure application working directory
WORKDIR /lco/noisegainreport

RUN apt-get update -y \
        && apt-get install --no-install-recommends -y less vim libpq-dev\
        && apt-get clean -y \
        && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python dependencies (libraries)
RUN pip --no-cache-dir install --upgrade pip \
        && pip install --upgrade --force-reinstall setuptools
        
# Install application code
COPY . .
RUN pip install .
