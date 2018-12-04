float random (vec2 st)
{
    return fract(sin(dot(st.xy, vec2(12.9898, 78.233)))* 43758.5453123);
}

float noise (vec2 st)
{
    vec2 i = floor(st);
    vec2 f = fract(st);

    // four corners in 2D of a tile
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));

    // smoothstep(0.,1.,f)
    vec2 u = f*f*(3.0-2.0*f);

    // mix 4 corners percentages
    return mix(a, b, u.x) +
            (c - a) * u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}

void main(void)
{
    vec2 st = gl_FragCoord.xy;
    float n = noise(st);
    gl_FragColor = vec4(vec3(n), 1.0);

    //gl_FragColor = vec4(0.8, 0.8, 0.0, 1.0);
}
