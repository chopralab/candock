#include <string.h>
#include <string.h>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>

#include <boost/filesystem.hpp>
#include "candock/drm/drm.hpp"

#include "candock/helper/logger.hpp"

using namespace std;
namespace candock {

#ifdef CANDOCK_DRM

static const unsigned int KEY_SIZE = 32;
static const unsigned int BLOCK_SIZE = 16;

/*
 * Encrypted file
 *
 *Unix Time
 *
 */

#include <openssl/conf.h>
#include <openssl/err.h>
#include <openssl/evp.h>

void handleErrors() {
    ERR_print_errors_fp(stderr);
    abort();
}

int encrypt(unsigned char* plaintext, int plaintext_len, unsigned char* key,
            unsigned char* iv, unsigned char* ciphertext) {
    EVP_CIPHER_CTX* ctx;

    int len;

    int ciphertext_len;

    /* Create and initialise the context */
    if (!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

    /* Initialise the encryption operation. IMPORTANT - ensure you use a key
     * and IV size appropriate for your cipher
     * In this example we are using 256 bit AES (i.e. a 256 bit key). The
     * IV size for *most* modes is the same as the block size. For AES this
     * is 128 bits */
    if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
        handleErrors();

    /* Provide the message to be encrypted, and obtain the encrypted output.
     * EVP_EncryptUpdate can be called multiple times if necessary
     */
    if (1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
        handleErrors();
    ciphertext_len = len;

    /* Finalise the encryption. Further ciphertext bytes may be written at
     * this stage.
     */
    if (1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len)) handleErrors();
    ciphertext_len += len;

    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);

    return ciphertext_len;
}

int decrypt(unsigned char* ciphertext, int ciphertext_len, unsigned char* key,
            unsigned char* iv, unsigned char* plaintext) {
    EVP_CIPHER_CTX* ctx;

    int len;

    int plaintext_len;

    /* Create and initialise the context */
    if (!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

    /* Initialise the decryption operation. IMPORTANT - ensure you use a key
     * and IV size appropriate for your cipher
     * In this example we are using 256 bit AES (i.e. a 256 bit key). The
     * IV size for *most* modes is the same as the block size. For AES this
     * is 128 bits */
    if (1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
        handleErrors();

    /* Provide the message to be decrypted, and obtain the plaintext output.
     * EVP_DecryptUpdate can be called multiple times if necessary
     */
    if (1 !=
        EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
        handleErrors();
    plaintext_len = len;

    /* Finalise the decryption. Further plaintext bytes may be written at
     * this stage.
     */
    if (1 != EVP_DecryptFinal_ex(ctx, plaintext + len, &len)) handleErrors();
    plaintext_len += len;

    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);

    return plaintext_len;
}

// Check DRM is designed to make sure the drm passes and if it does then
// create a new file with a later time stamp to prevent further tampering
// Passed DRM if returns true
// Failed DRM if returns false
bool drm::check_drm(const std::string& key_location) {
    bool passed_drm = true;

    /* Initialise the library */
    ERR_load_crypto_strings();
    OpenSSL_add_all_algorithms();
    OPENSSL_config(NULL);

    /* A 256 bit key */
    unsigned char* key = (unsigned char*)"25280559686690074811771763602098";

    /* A 128 bit IV */
    unsigned char* iv = (unsigned char*)"1226610434284021";

    /*
     *    Check to see if the encrypted file exists.
     *    If it does decrypt it.
     */

    if (boost::filesystem::exists(key_location)) {
        string readEncrypted;
        ifstream readEncryptedFile(key_location);
        getline(readEncryptedFile, readEncrypted);
        readEncryptedFile.close();

        /* Message to be encrypted */
        unsigned char encrypted[128];
        strcpy((char*)encrypted, readEncrypted.c_str());

        /* Buffer for decryptedtext. Ensure the buffer is long enough for the
         * decryptedtext which may be longer than the plaintext, dependant on
         * the
         * algorithm and mode
         */
        unsigned char decryptedtext[128];

        /* Decrypt the ciphertext */
        decrypt(encrypted, 16, key, iv, decryptedtext);

        /* Do something useful with the ciphertext here */
        // printf("Ciphertext is:\n");
        // BIO_dump_fp (stdout, (const char *)decryptedtext, decryptedtext_len);

        /* Add a NULL terminator. We are expecting printable text */
        // decryptedtext[decryptedtext_len] = '\0';

        // cout << "Number " << decryptedtext << endl;
        // unsigned int num =  static_cast<unsigned int>(decryptedtext);

        const char* end_date = "1504936752";
        const char* start_date = "1494809088";

        if (strncmp((const char*)decryptedtext, start_date, 10) < 0 ||
            strncmp((const char*)decryptedtext, end_date, 10) > 0) {
            passed_drm = false;
        }

    } else {
        // This is what happens if the file does not exist
        log_error << "Error: No Key Found!"
                  << "\n";
        passed_drm = false;
    }

    string line = std::to_string(static_cast<uintmax_t>(std::time(NULL)));

    /* Message to be encrypted */
    unsigned char plaintext[128];

    /* Buffer for ciphertext. Ensure the buffer is long enough for the
     * ciphertext which may be longer than the plaintext, dependant on the
     * algorithm and mode
     */
    unsigned char ciphertext[128];
    strcpy((char*)plaintext, line.c_str());
    /* Encrypt the plaintext */
    encrypt(plaintext, strlen((char*)plaintext), key, iv, ciphertext);

    if (passed_drm) {
        ofstream myfile3;
        myfile3.open(key_location);
        myfile3 << ciphertext;
        myfile3.close();
    }

    /* Clean up */
    EVP_cleanup();
    ERR_free_strings();

    return passed_drm;
}

#else

bool drm::check_drm(const std::string&) { return true; }

#endif /*CANDOCK_DRM*/
}
