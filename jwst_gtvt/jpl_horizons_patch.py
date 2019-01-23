"""Patch for vectors class, necessary until astropy PR is merged (hopefully.)
"""

from collections import OrderedDict
from numpy import nan as nan
from numpy import isnan
from numpy import ndarray

from astropy.time import Time

from astroquery.jplhorizons import HorizonsClass
from astroquery.jplhorizons import conf

def vectors_async_full(self, get_query_payload=False, refplane='ecliptic',
                    closest_apparition=False, no_fragments=False,
                    get_raw_response=False, cache=True):
    
    URL = conf.horizons_server

    # check for required information
    if self.id is None:
        raise ValueError("'id' parameter not set. Query aborted.")
    if self.location is None:
        self.location = '500@10'
    if self.epochs is None:
        self.epochs = Time.now().jd

    # assemble commandline based on self.id_type
    commandline = str(self.id)

    if self.id_type in ['designation', 'name',
                        'asteroid_name', 'comet_name']:
        commandline = ({'designation': 'DES=',
                        'name': 'NAME=',
                        'asteroid_name': 'ASTNAM=',
                        'comet_name': 'COMNAM='}[self.id_type] +
                        commandline)
    if self.id_type in ['smallbody', 'asteroid_name',
                        'comet_name', 'designation']:
        commandline += ';'
        if isinstance(closest_apparition, bool):
            if closest_apparition:
                commandline += ' CAP;'
        else:
            commandline += ' CAP{:s};'.format(closest_apparition)
        if no_fragments:
            commandline += ' NOFRAG;'

    if isinstance(self.location, dict):
        raise ValueError(('cannot use topographic position in state'
                            'vectors query'))

    # configure request_payload for ephemerides query
    request_payload = OrderedDict([
        ('batch', 1),
        ('TABLE_TYPE', 'VECTORS'),
        ('OUT_UNITS', 'AU-D'),
        ('COMMAND', '"' + commandline + '"'),
        ('CENTER', ("'" + str(self.location) + "'")),
        ('CSV_FORMAT', ('"YES"')),
        ('REF_PLANE', {'ecliptic': 'ECLIPTIC', 'earth': 'FRAME',
                        'body': "'BODY EQUATOR'"}[refplane]),
        ('REF_SYSTEM', 'J2000'),
        ('TP_TYPE', 'ABSOLUTE'),
        ('LABELS', 'YES'),
        ('OBJ_DATA', 'YES')]
    )

    # parse self.epochs
    if isinstance(self.epochs, (list, tuple, ndarray)):
        request_payload['TLIST'] = "\n".join([str(epoch) for epoch in
                                                self.epochs])
    elif type(self.epochs) is dict:
        if ('start' not in self.epochs or 'stop' not in self.epochs or
                'step' not in self.epochs):
            raise ValueError("'epochs' must contain start, " +
                                "stop, step")
        request_payload['START_TIME'] = self.epochs['start']
        request_payload['STOP_TIME'] = self.epochs['stop']
        request_payload['STEP_SIZE'] = self.epochs['step']

    else:
        # treat epochs as a list
        request_payload['TLIST'] = str(self.epochs)

    self.query_type = 'vectors'

    # return request_payload if desired
    if get_query_payload:
        return request_payload

    # set return_raw flag, if raw response desired
    if get_raw_response:
        self.return_raw = True

    # query and parse
    response = self._request('GET', URL, params=request_payload,
                                timeout=self.TIMEOUT, cache=cache)
    self.uri = response.url

    return response