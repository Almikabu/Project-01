from flask import Flask, request, render_template, send_file
import DNAencrypt
import DNAdecrypt
from utils import str2bin, bin2str
import random
import string
import io

app = Flask(__name__)

# Global variable to store the most recent decryption key
global_decryption_key = None

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/encrypt', methods=['POST'])
def encrypt():
    global global_decryption_key
    text_to_encrypt = request.form['text_to_encrypt']
    if not text_to_encrypt:
        return render_template('index.html', error='Please enter text to encrypt')
    
    try:
        DNAencrypt.set_globals()
        key = str2bin(''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(16)))
        DNAencrypt.generate_pre_processing_tables()
        DNAencrypt.generate_mutation_tables()
        
        encrypted_text = DNAencrypt.dnaEncrypt(text_to_encrypt, key)
        global_decryption_key = DNAencrypt.decryption_key
        
        return render_template('index.html', encrypted_text=encrypted_text)
    except Exception as e:
        return render_template('index.html', error=f'Encryption error: {str(e)}')

@app.route('/download_key')
def download_key():
    global global_decryption_key
    if global_decryption_key is None:
        return "No key available. Please encrypt a message first.", 400
    
    return send_file(
        io.BytesIO(global_decryption_key.encode()),
        mimetype='text/plain',
        as_attachment=True,
        download_name='decryption_key.txt'
    )

@app.route('/decrypt', methods=['POST'])
def decrypt():
    text_to_decrypt = request.form['text_to_decrypt']
    if 'key_file' not in request.files:
        return render_template('index.html', error='No key file uploaded')
    key_file = request.files['key_file']
    if key_file.filename == '':
        return render_template('index.html', error='No key file selected')
    
    try:
        decryption_key = key_file.read().decode('utf-8')
        decrypted_text = DNAdecrypt.dnaDecrypt(text_to_decrypt, decryption_key)
        return render_template('index.html', decrypted_text=decrypted_text)
    except Exception as e:
        return render_template('index.html', error=f'Decryption error: {str(e)}')

if __name__ == '__main__':
    app.run(debug=True)