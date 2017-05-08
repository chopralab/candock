#include <string.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>

#include <boost/filesystem.hpp>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include "drm.hpp"

using namespace std;
static const unsigned int KEY_SIZE = 32;
static const unsigned int BLOCK_SIZE = 16;

/*
 *Encrypted file
 *
 *Unix Time
 *
 */

 void handleErrors()
{
  ERR_print_errors_fp(stderr);
  abort();
}


 int decrypt(unsigned char *ciphertext, int ciphertext_len, unsigned char *key, unsigned char *iv, unsigned char *plaintext)
{
  EVP_CIPHER_CTX *ctx;

  int len;

  int plaintext_len;

  /* Create and initialise the context */
  if(!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

  /* Initialise the decryption operation. IMPORTANT - ensure you use a key
   * and IV size appropriate for your cipher
   * In this example we are using 256 bit AES (i.e. a 256 bit key). The
   * IV size for *most* modes is the same as the block size. For AES this
   * is 128 bits */
  if(1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    handleErrors();

  /* Provide the message to be decrypted, and obtain the plaintext output.
   * EVP_DecryptUpdate can be called multiple times if necessary
   */
  if(1 != EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
    handleErrors();
  plaintext_len = len;

  /* Finalise the decryption. Further plaintext bytes may be written at
   * this stage.
   */
  if(1 != EVP_DecryptFinal_ex(ctx, plaintext + len, &len)) handleErrors();
  plaintext_len += len;

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);

  return plaintext_len;
}



//Check DRM is designed to make sure the drm passes and if it does then
//create a new file with a later time stamp to prevent further tampering
//Passed DRM if returns true
//Failed DRM if returns false
 bool drm::check_drm() {

  /* Initialise the library */
  ERR_load_crypto_strings();
  OpenSSL_add_all_algorithms();
  OPENSSL_config(NULL);

   /* A 256 bit key */
  unsigned char *key = (unsigned char *)"01234567890123456789012345678901";

  /* A 128 bit IV */
  unsigned char *iv = (unsigned char *)"0123456789012345";

  /*
    Check to see if the encrypted file exists.
    If it does decrypt it.
  */
  if(boost::filesystem::exists("/home/brandon_stewart/.candock")) {

    string readEncrypted;
    ifstream readEncryptedFile ("/home/brandon_stewart/.candock");
    getline (readEncryptedFile,readEncrypted);
    cout << readEncrypted << '\n';
    readEncryptedFile.close();
    
    /* Message to be encrypted */
    unsigned char encrypted[128];
    strcpy( (char*) encrypted, readEncrypted.c_str()); 

    /* Buffer for decryptedtext. Ensure the buffer is long enough for the
     * decryptedtext which may be longer than the plaintext, dependant on the
     * algorithm and mode
     */
    unsigned char decryptedtext[128];


    /* Decrypt the ciphertext */
    int decryptedtext_len = decrypt(encrypted, 16, key, iv,
    decryptedtext);

    if((long long int) decryptedtext > 1504936752 || (long long int) std::time(NULL) < 1494226769) {
        /* Clean up */
        EVP_cleanup();
        ERR_free_strings();
        return false;
    }

/* Do something useful with the ciphertext here */
  printf("Ciphertext is:\n");
  BIO_dump_fp (stdout, (const char *)decryptedtext, decryptedtext_len);


    /* Add a NULL terminator. We are expecting printable text */
    decryptedtext[decryptedtext_len] = '\0';

    cout << "Number " << decryptedtext << endl;


  }
  else {
      cout << "Hello" << endl;
     if(std::time(NULL) < 1494226769 || std::time(NULL) > 1504936752) {
        /* Clean up */
        EVP_cleanup();
        ERR_free_strings();
        return false;
     }
  }

  /* Buffer for the decrypted text */
 /* unsigned char ciphertext[128];
  unsigned char plaintext[128];
  
  plaintext = static_cast<unsigned char*> std::time(NULL);

*/
  /* Encrypt the plaintext */
 /* int ciphertext_len = encrypt (plaintext, strlen ((char *)plaintext), key, iv,
                            ciphertext);

  ofstream writeEncrypted;
  writeEncrypted.open ("~/.candock");
  writeEncrypted << ciphertext;
  writeEncrypted.close();
*/

  /* Clean up */
  EVP_cleanup();
  ERR_free_strings();

  return true;
}










 int encrypt(unsigned char *plaintext, int plaintext_len, unsigned char *key,
  unsigned char *iv, unsigned char *ciphertext)
{
  EVP_CIPHER_CTX *ctx;

  int len;

  int ciphertext_len;

  /* Create and initialise the context */
  if(!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

  /* Initialise the encryption operation. IMPORTANT - ensure you use a key
   * and IV size appropriate for your cipher
   * In this example we are using 256 bit AES (i.e. a 256 bit key). The
   * IV size for *most* modes is the same as the block size. For AES this
   * is 128 bits */
  if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    handleErrors();

  /* Provide the message to be encrypted, and obtain the encrypted output.
   * EVP_EncryptUpdate can be called multiple times if necessary
   */
  if(1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
    handleErrors();
  ciphertext_len = len;

  /* Finalise the encryption. Further ciphertext bytes may be written at
   * this stage.
   */
  if(1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len)) handleErrors();
  ciphertext_len += len;

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);

  return ciphertext_len;
}

