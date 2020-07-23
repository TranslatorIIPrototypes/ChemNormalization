# Chemical Normalization

Chemical normalization takes a chemical substance (CURIE) and returns a list of similar chemicals.

Each chemical substance is defined by a SMILES value. In this application we simplify each chemical substance''s SMILES using RDKit and then group these chemicals by their SMILES values.
    
## Deployment

Assume a Redis database running on localhost port 6379.

### Local

```bash
pip install -r requirements.txt
export REDIS_HOST=localhost
export REDIS_PORT=6379
./main.py --port 6380
```

### Docker

```bash
docker build -t redis_rest .
docker run \
    -p 6380:6380 \
    --env REDIS_HOST=host.docker.internal \
    --env REDIS_PORT=6379 \
    redis_rest --port 6380
```

### Docker compose

To run both Redis database and REST interface, a docker-compose file is provided:

```bash
docker-compose up
```

## usage

    <https://chemnormalization.renci.org/apidocs/>