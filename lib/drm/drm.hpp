#ifndef DRM_H
#define DRM_H

class drm{

    private:
        static void encrypt();
        static bool decrypt();

    public:
        static bool check_drm();

};

#endif
