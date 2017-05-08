#include <iostream>
#include <string>
#include <memory>
#include <limits>
#include <stdexcept>

#include <openssl/evp.h>
#include <openssl/rand.h>



static const unsigned int KEY_SIZE = 32;
static const unsigned int BLOCK_SIZE = 16;

//Passed DRM if returns true
//Failed DRM if returns false
static bool drm::check_drm() {
    if(decrpyt()) {
        return false;
    }

    return true;
}

static void drm::encrypt() {

}

static bool drm::decrypt() {

}
