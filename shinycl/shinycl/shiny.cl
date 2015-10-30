kernel void init(global float *t, global float *t_next, uint nx, uint ny) {
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);
    
    if (i < nx - 1 && j < ny - 1) {
       size_t k_ij   = i*ny + j;
    
       t[k_ij] = 0;
       t_next[k_ij] = 0;
    }
}

float position(size_t index, float slope, float offset) {
    return index*slope + offset;
}

float2 get_b(float x, float y) {
    float2 b;
    b.x =  M_PI_F*cospi(x)*sinpi(y);
    b.y = -M_PI_F*cospi(y)*sinpi(x);
    
    float bmod = sqrt(b.x*b.x + b.y*b.y);
    
    return b/bmod;
}

float2 dDotGradT(float kper, float kpar, float2 b, float2 gT) {
    float dxx = b.x*b.x*(kpar - kper) + kper;
    float dxy = b.x*b.y*(kpar - kper);
    float dyy = b.y*b.y*(kpar - kper) + kper;
    
    float2 q;
    q.x = dxx*gT.x + dxy*gT.y;
    q.y = dxy*gT.x + dyy*gT.y;
    
    return q;
}

kernel void stepTimeGunterSym(global float *t, global float *t_next,
                              uint nx, uint ny,
                              float xmin, float ymin,
                              float dt, float dx, float dy,
                              float kper, float kpar) {
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);
    
    if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
        float x = position(i, dx, xmin); // Check the blance of data access to flops this can be preprocessed.
        float y = position(j, dy, ymin); // Check the blance of data access to flops this can be preprocessed.
        
        float xp = x + dx/2; // Check the blance of data access to flops this can be preprocessed.
        float xm = x - dx/2; // Check the blance of data access to flops this can be preprocessed.
        float yp = y + dy/2; // Check the blance of data access to flops this can be preprocessed.
        float ym = y - dy/2; // Check the blance of data access to flops this can be preprocessed.
        
        size_t k_ipjp = (i + 1)*ny + (j + 1); // Check the blance of data access to flops this can be preprocessed.
        size_t k_ipj  = (i + 1)*ny + j;       // Check the blance of data access to flops this can be preprocessed.
        size_t k_ipjm = (i + 1)*ny + (j - 1); // Check the blance of data access to flops this can be preprocessed.
        size_t k_ijp  = i*ny + (j + 1);       // Check the blance of data access to flops this can be preprocessed.
        size_t k_ij   = i*ny + j;             // Check the blance of data access to flops this can be preprocessed.
        size_t k_ijm  = i*ny + (j - 1);       // Check the blance of data access to flops this can be preprocessed.
        size_t k_imjp = (i - 1)*ny + (j + 1); // Check the blance of data access to flops this can be preprocessed.
        size_t k_imj  = (i - 1)*ny + j;       // Check the blance of data access to flops this can be preprocessed.
        size_t k_imjm = (i - 1)*ny + (j - 1); // Check the blance of data access to flops this can be preprocessed.
        
        float2 gT_ipjp;
        gT_ipjp.x = (t[k_ipjp] + t[k_ipj]  - t[k_ijp]  - t[k_ij])/(2*dx);
        gT_ipjp.y = (t[k_ijp]  + t[k_ipjp] - t[k_ipj]  - t[k_ij])/(2*dy);
        
        float2 gT_ipjm;
        gT_ipjm.x = (t[k_ipj]  + t[k_ipjm] - t[k_ij]   - t[k_ijm])/(2*dx);
        gT_ipjm.y = (t[k_ij]   + t[k_ipj]  - t[k_ipjm] - t[k_ijm])/(2*dy);
        
        float2 gT_imjp;
        gT_imjp.x = (t[k_ijp]  + t[k_ij]   - t[k_imjp] - t[k_imj])/(2*dx);
        gT_imjp.y = (t[k_imjp] + t[k_ijp]  - t[k_ij]   - t[k_imj])/(2*dy);
        
        float2 gT_imjm;
        gT_imjm.x = (t[k_ij]   + t[k_ijm]  - t[k_imj]  - t[k_imjm])/(2*dx);
        gT_imjm.y = (t[k_imj]  + t[k_ij]   - t[k_ijm]  - t[k_imjm])/(2*dy);
        
        float2 q_ipjp = -dDotGradT(kper, kpar, get_b(xp, yp), gT_ipjp); // Check the blance of data access to flops b can be preprocessed.
        float2 q_ipjm = -dDotGradT(kper, kpar, get_b(xp, ym), gT_ipjm); // Check the blance of data access to flops b can be preprocessed.
        float2 q_imjp = -dDotGradT(kper, kpar, get_b(xm, yp), gT_imjp); // Check the blance of data access to flops b can be preprocessed.
        float2 q_imjm = -dDotGradT(kper, kpar, get_b(xm, ym), gT_imjm); // Check the blance of data access to flops b can be preprocessed.
        
        float divq = (q_ipjp.x + q_ipjm.x - q_imjp.x - q_imjm.x)/(2*dx)
                   + (q_ipjp.y + q_imjp.y - q_imjm.y - q_ipjm.y)/(2*dy);
        
        float q = 2*M_PI_F*M_PI_F*cospi(x)*cospi(y); // Check the blance of data access to flops this can be preprocessed.
        
        t_next[k_ij] = dt*(-divq + q);
    }
}

kernel void update(global float *t, global float *t_next, uint nx, uint ny) {
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);
    
    if (i < nx - 1 && j < ny - 1) {
       size_t k_ij = i*ny + j;
    
       t[k_ij] += t_next[k_ij];
    }
}

kernel void analytic(global float *t,
                     uint nx, uint ny,
                     float xmin, float ymin,
                     float time, float dx, float dy,
                     float kper) {
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);
    
    if (i < nx - 1 && j < ny - 1) {
        size_t k_ij = i*ny + j;
        
        float x = position(i, dx, xmin);
        float y = position(j, dy, ymin);
        
        float psi = cospi(x)*cospi(y);
        
        t[k_ij] = (1.0 - exp(-2*kper*M_PI_F*M_PI_F*time))/kper*psi;
    }
}
