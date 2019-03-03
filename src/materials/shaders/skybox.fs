//precision mediump float;
//precision mediump int;

uniform vec3 sunDirection;

#define atmosphereRadius  6420e3
#define earthRadius  6360e3
#define Hr  7994.0
#define Hm  1200.0
#define betaR vec3(3.8e-6, 13.5e-6, 33.1e-6)
#define betaM vec3(21e-6, 21e-6, 21e-6)

varying vec3 pos;

vec3 solveQuadratic(float a, float b, float c, float x1, float x2)
{
    if (b == 0.0) {
        // Handle special case where the the two vector ray.dir and V are perpendicular
        // with V = ray.orig - sphere.centre
        if (a == 0.0) return vec3(0.0, x1, x2);
        x1 = 0.0; x2 = sqrt(-c / a);
        return vec3(1.0, x1, x2);
    }
    float discr = b * b - 4.0 * a * c;

    if (discr < 0.0) return vec3(0.0, x1, x2);

    float q = (b < 0.0) ? -0.5 * (b - sqrt(discr)) : -0.5 * (b + sqrt(discr));
    x1 = q / a;
    x2 = c / q;

    return vec3(1.0, x1, x2);
}


vec3 raySphereIntersect(vec3 orig, vec3 dir, float radius, float t0, float t1)
{
    // They ray dir is normalized so A = 1 
    float A = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;
    float B = 2.0 * (dir.x * orig.x + dir.y * orig.y + dir.z * orig.z);
    float C = orig.x * orig.x + orig.y * orig.y + orig.z * orig.z - radius * radius;

    vec3 r = solveQuadratic(A, B, C, t0, t1);
    if (r.x < 1.0) return vec3(0.0, r.y, r.z); // false

    if (r.y > r.z) return vec3(1.0, r.z, r.y);

    return vec3(1.0, r.y, r.z); // true
}

#define M_PI (3.14159265358979323846) 


vec3 computeIncidentLight(vec3 orig, vec3 dir, vec3 sunDirection, float tmin, float tmax) {
    float t0, t1;
    vec3 r = raySphereIntersect(orig, dir, atmosphereRadius, t0, t1);
    t0 = r.y;
    t1 = r.z;
    if (r.x < 1.0 || t1 < 0.0) return vec3(1.0, 0.0, 0.0);
    if (t0 > tmin && t0 > 0.0) tmin = t0;
    if (t1 < tmax) tmax = t1;
    float numSamples = 16.0;
    float numSamplesLight = 8.0;
    float segmentLength = (tmax - tmin) / numSamples;
    float tCurrent = tmin;
    vec3 sumR = vec3(0.0, 0.0, 0.0);
    vec3 sumM = vec3(0.0, 0.0, 0.0); // mie and rayleigh contribution
    float opticalDepthR = 0.0;
    float opticalDepthM = 0.0;
    float mu = dot(dir, sunDirection); // mu in the paper which is the cosine of the angle between the sun direction and the ray direction
    float phaseR = 3.0 / (16.0 * M_PI) * (1.0 + mu * mu);
    float g = 0.76;
    float phaseM = 3.0 / (8.0 * M_PI) * ((1.0 - g * g) * (1.0 + mu * mu)) / ((2.0 + g * g) * pow(1.0 + g * g - 2.0 * g * mu, 1.5));
    for (int i = 0; i < 16; i++) {
        vec3 samplePosition = orig + (tCurrent + segmentLength * 0.5) * dir;
        float height = length(samplePosition) - earthRadius;
        // compute optical depth for light
        float hr = exp(-height / Hr) * segmentLength;
        float hm = exp(-height / Hm) * segmentLength;
        opticalDepthR += hr;
        opticalDepthM += hm;
        // light optical depth
        float t0Light, t1Light;
        vec3 r = raySphereIntersect(samplePosition, sunDirection, atmosphereRadius, t0Light, t1Light);
        t0Light = r.y;
        t1Light = r.z;
        float segmentLengthLight = t1Light / numSamplesLight, tCurrentLight = 0.0;
        float opticalDepthLightR = 0.0;
        float opticalDepthLightM = 0.0;
        //int j;
        bool done = true;
        for (int j = 0; j < 8; j++) {
            vec3 samplePositionLight = samplePosition + (tCurrentLight + segmentLengthLight * 0.5) * sunDirection;
            float heightLight = length(samplePositionLight) - earthRadius;
            if (heightLight < 0.0) {
                done = false;
                break;
            }
            opticalDepthLightR += exp(-heightLight / Hr) * segmentLengthLight;
            opticalDepthLightM += exp(-heightLight / Hm) * segmentLengthLight;
            tCurrentLight += segmentLengthLight;
        }
        //if (j == 8) {
        if(done){
            vec3 tau = betaR * (opticalDepthR + opticalDepthLightR) + betaM * 1.1 * (opticalDepthM + opticalDepthLightM);
            vec3 attenuation = vec3(exp(-tau.x), exp(-tau.y), exp(-tau.z));
            sumR += attenuation * hr;
            sumM += attenuation * hm;
        }
        //}
        tCurrent += segmentLength;
    }

    // We use a magic number here for the intensity of the sun (20). We will make it more
    // scientific in a future revision of this lesson/code
    return (sumR * betaR * phaseR + sumM * betaM * phaseM) * 20.0;
} 
float infinity = 1.0 / 0.0;
float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

void main(){
    vec3 orig = vec3(0.0, earthRadius + 1.0, 0.0);
    vec3 SunPosition = sunDirection;
	vec3 sunVector = vec3(
		-sin(SunPosition.x)*cos(SunPosition.y),
		sin(SunPosition.y),
		cos(SunPosition.x)*cos(SunPosition.y)
	);
	vec3 c = computeIncidentLight(orig, normalize(pos), sunVector, 0.0, infinity);
	gl_FragColor = vec4(c, 1.0);
}
