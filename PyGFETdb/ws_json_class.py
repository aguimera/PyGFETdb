#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
#
import json
import traceback
import requests
import base64

WS_ENTRY = "http://opter6.cnm.es/ws/pyfet/pyfet.php"


class ws_json:
    '''
        Paràmetres:
            query: Les dades en json que enviarem al servei web
        Retorn: Array amb el resultat:
            'error' => [true|false]
            'errorCode' => Codi de l'error (0) si no hi ha error
            'message' => En cas d'error, el missatge d'error.
                Si no hi ha error, aquest camp no existeix
                'result' => Resultat de la crida
    '''
    @staticmethod
    def send(query):
        data_string = json.dumps(query)
        headers = {'Content-Type': 'application/json'}
        response = requests.post(WS_ENTRY,
                                 headers=headers,
                                 data=data_string)
        try:
            result = json.loads(response.text)
        except:
            result = {'error': 1,
                      'errorCode': 1,
                      'message': 'Cannot exec json request'}
            return result

        try:
            error = result['error']
        except:
            result = {'error': 1,
                      'errorCode': 2,
                      'message': 'malformed answer ...'}

        try:
            result['result'] = ws_json.ws_unpack(result['result'])
        except:
            traceback.print_exc()

        return(result)

    '''
        Aquesta funció composa l'array amb els paràmetres que enviarem al servei web
        Paràmetres:
            module: El mòdul al que fem referència. Pot set pyfet o pyegnite
            data: La crida sql que volem executar
        Retorn: El resultat de la crida al servei web
            Array amb el resultat:
            'error' => [true|false]
            'errorCode' => Codi de l'error (0) si no hi ha error
            'message' => En cas d'error, el missatge d'error.
                        Si no hi ha error, aquest camp no existeix
            'result' => Resultat de la crida
    '''
    @staticmethod
    def call(module, data):
        username = 'py-imb-dbtools'
        password = '6345HGdfkfmjk'
        query = {'username': username,
                 'password': password,
                 'module': module,
                 'oper': 'query',
                 'sql': data}
        response = ws_json.send(query)
        return(response)

    '''
        Aquesta funció retorna una llista on cada element es un array (dictionary)
        Si un valor en concret el volem en binary, podem fer bytesarray(valor)
        Si un valor en concret el volem en string, podem fer valor.decode()
        (valor és cada un dels elements de row[j]) ...
        Paràmetres:
            data => llista d'entrada on cada element és un array associatiu amb el valor en b64
            Retorn: La llista amb els valors decodificats. 
                    Cada element de l'array que forma cada element de la llista
                    es retorna en 'bytes'
	'''
    @staticmethod
    def ws_unpack(data):
        retVal = []
        for i in data:
            row = {}
            for k, v in i.items():
                row[k] = base64.b64decode(v)
            retVal.append(row)
        return(retVal)
