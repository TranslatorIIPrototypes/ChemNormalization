"""Sanic server for Redis-REST with referencing."""
import json
import os

import aioredis
from sanic import Sanic, response

from r3.apidocs import bp as apidocs_blueprint

app = Sanic()
app.config.ACCESS_LOG = False
app.blueprint(apidocs_blueprint)

redis_host = os.environ.get('REDIS_HOST', 'localhost')
redis_port = os.environ.get('REDIS_PORT', 6379)


async def startup_connections(app, loop):
    """Start up Redis connection."""
    app.redis_connection0 = await aioredis.create_redis_pool(
        f'redis://{redis_host}:{redis_port}', db=0)
    app.redis_connection1 = await aioredis.create_redis_pool(
        f'redis://{redis_host}:{redis_port}', db=1)


async def shutdown_connections(app, loop):
    """Shut down Redis connection."""
    app.redis_connection0.close()
    await app.redis_connection0.wait_closed()
    app.redis_connection1.close()
    await app.redis_connection1.wait_closed()


@app.route('/get_similar_chemicals')
async def get_handler(request):
    """Get value(s) for key(s).

    Use GET for single key, MGET for multiple.
    """
    if isinstance(request.args['curie'], list):
        references = await app.redis_connection0.mget(*request.args['curie'], encoding='utf-8')
        references_nonnan = [reference for reference in references if reference is not None]
        if not references_nonnan:
            return response.text('No matches found for specified CURIE(s).', status=404)
        values = await app.redis_connection1.mget(*references_nonnan, encoding='utf-8')
        values = [json.loads(value) if value is not None else None for value in values]
        dereference = dict(zip(references_nonnan, values))
        return response.json({
            key: dereference[reference] if reference is not None else None
            for key, reference in zip(request.args['curie'], references)
        })
    else:
        reference = await app.redis_connection0.get(request.args['curie'], encoding='utf-8')
        if reference is None:
            return response.json({request.args['curie']: None})
        value = await app.redis_connection1.get(reference, encoding='utf-8')
        value = json.loads(value) if value is not None else None
        return response.json({request.args['curie']: value})


app.register_listener(startup_connections, 'after_server_start')
app.register_listener(shutdown_connections, 'before_server_stop')
