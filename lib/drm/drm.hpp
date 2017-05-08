#ifndef DRM_H
#define DRM_H

class drm {

   
        static int encrypt(unsigned char *plaintext, int plaintext_len, unsigned char *key, unsigned char *iv, unsigned char *ciphertext);
        //static int decrypt(unsigned char *ciphertext, int ciphertext_len, unsigned char *key, unsigned char *iv, unsigned char *plaintext);
       static  void handleErrors();

    public:
        static bool check_drm();

};

#endif
