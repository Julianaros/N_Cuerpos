# SIMULACIÓN DEL PROBLEMA DE LOS N-CUERPOS COMPARANDO TÉCNICAS DE PARALELIZACIÓN

Este proyecto se busca simular el problema de los N-cuerpos, un problema fundamental en física computacional para simulaciones con aplicación en astronomía y dinámica molecular por mencionar algunos. A medida que la cantidad de cuerpos aumenta el problema se complejiza a tal grado de necesitar paralelización para obtener resultados en tiempos razonables. 

###  Ecuaciones Fundamentales

Partimos de la **ecuación de Newton**, donde cada partícula experimenta la fuerza

\( F_i = \sum_{j=1,\, j \ne i}^{N} G \frac{m_i m_j}{|\vec{r}_{ij}|^3} \vec{r}_{ij} \)

y usando la **segunda Ley de Newton** es posible conocer la aceleración de las partículas

\( \vec{a}_i = \frac{\vec{F}_i}{m_i} = \sum_{\substack{j=1 \\ j \ne i}}^{N} G \frac{m_j}{|\vec{r}_{ij}|^3} \vec{r}_{ij} \)

donde
- \( \vec{r}_{ij} = \vec{r}_j - \vec{r}_i \) es el vector de distancia  
- \( G \) es la constante gravitacional (normalizada a 1)  
- \( m_i, m_j \) son las masas de las partículas  
