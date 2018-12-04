void main(void)
{
    vec4 osVert = gl_Vertex;
    vec4 esVert = gl_ModelViewMatrix * osVert;
    vec4 csVert = gl_ProjectionMatrix * esVert;
    gl_Position = csVert;

    //gl_Position = gl_ModelViewProjectionMatrix * glVertex;
    gl_TexCoord[0] = gl_MultiTexCoord0;
}
