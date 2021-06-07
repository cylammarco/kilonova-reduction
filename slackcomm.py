import json
import requests
from encode import decode

# post message
def post_message_to_slack(webhook_url, message, channels="test"):

    response = requests.post(url=webhook_url,
                             data=json.dumps({'text': message}),
                             headers={'Content-Type': 'application/json'})

    return response


# post file
def post_file_to_slack(slack_token,
                       filename,
                       filetype="auto",
                       channels="test"):

    response = requests.post(url='https://slack.com/api/files.upload',
                             data={
                                 'token': slack_token,
                                 'filename': filename,
                                 'filetype': filetype,
                                 'channels': channels
                             },
                             files={'file': (filename, open(filename, 'rb'))})

    return response
