import base64

def encode(key, string):

    encoded_chars = []

    for i in range(len(string)):

        key_c = key[i % len(key)]
        encoded_c = chr(ord(string[i]) + ord(key_c) % 128)
        encoded_chars.append(encoded_c)

    encoded_string = "".join(encoded_chars)
    arr2 = bytes(encoded_string, 'utf-8')

    return base64.urlsafe_b64encode(arr2)

def decode(key, string):

    encoded_chars = []
    string = base64.urlsafe_b64decode(string)
    string = string.decode('utf-8')

    for i in range(len(string)):

        key_c = key[i % len(key)]
        encoded_c = chr(ord(string[i]) - ord(key_c) % 128)
        encoded_chars.append(encoded_c)

    encoded_string = "".join(encoded_chars)

    return encoded_string
