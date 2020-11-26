import base64
import os
from cryptography.fernet import Fernet

def encode(password, string):

    if len(password) > 32:
        raise ValueError('Password has to be shorter than 32 charaters.')
    len_diff = 32 - len(password)
    password += "0" * len_diff

    key = base64.urlsafe_b64encode(bytes(password, 'utf-8'))
    Fkey = Fernet(key)
    token = Fkey.encrypt(bytes(string, 'utf-8'))

    return token

def decode(password, token):

    if len(password) > 32:
        raise ValueError('Password has to be shorter than 32 charaters.')
    len_diff = 32 - len(password)
    password += "0" * len_diff
    
    key = base64.urlsafe_b64encode(bytes(password, 'utf-8'))
    Fkey = Fernet(key)
    message = Fkey.decrypt(bytes(token, 'utf-8'))

    return message

