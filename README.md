# SIMULACIÓN DEL PROBLEMA DE LOS N-CUERPOS COMPARANDO TÉCNICAS DE PARALELIZACIÓN

Este proyecto se busca simular el problema de los N-cuerpos, un problema fundamental en física computacional para simulaciones con aplicación en astronomía y dinámica molecular por mencionar algunos. A medida que la cantidad de cuerpos aumenta el problema se complejiza a tal grado de necesitar paralelización para obtener resultados en tiempos razonables. 

###  Ecuaciones Fundamentales

Partimos de la **ecuación de Newton**, donde cada partícula experimenta la fuerza

![image](https://github.com/user-attachments/assets/39066bb2-dd71-47d6-999e-737d482604ab)

y usando la **segunda Ley de Newton** es posible conocer la aceleración de las partículas

![image](https://github.com/user-attachments/assets/405375c4-d0eb-4aaf-a289-cb59da2112f4)

donde

![image](https://github.com/user-attachments/assets/80c2942d-02c0-4e04-a727-178b3b1e2715)

### Solución por el algoritmo de Verlet

Este conocido método matemático permite calcular la evolución temporal de la posición de una partícula bajo la acción de una fuerza conocida usando posiciones pasadas previamente definidas o calculadas. 

![image](https://github.com/user-attachments/assets/2f2e0118-619e-4cf8-bcc0-b52cc0f218ae)

### Paralelización implementada

En este caso, para la comparación de eficiencias y rendimiento, se evalúan los siguientes paradigmas de paralelización:
- OpenMP
- Pthreads
- MPI

### Referencias
1. Aversa, R., Di Martino, B., Mazzocca, N., \& Venticinque, S. (2005). Performance analysis of hybrid OpenMP/MPI N-body application. In B. M. Chapman (Ed.), Shared memory parallel programming with Open MP (WOMPAT 2004) (Lecture Notes in Computer Science, Vol. 3349, pp. 11–22). Springer. https://doi.org/10.1007/978-3-540-31832-3\_2
2. Kang, S., \& Lee, S. (2015). Performance comparison of OpenMP, MPI, and MapReduce in practical problems. Advances in Multimedia, 2015, 1–9. https://doi.org/10.1155/2015/575687
3. Kuzmin, N., Sirotin, D., \& Khoperskov, A. (2024). Efficiency of parallel computations of gravitational forces by treecode method in N-body models. Mathematical Physics and Computer Simulation, 27(4), 39–55. https://doi.org/10.15688/mpcm.jvolsu.2024.4.4
4. Truong Vinh Truong Duy, T., Yamazaki, K., Ikegami, K., \& Oyanagi, S. (2012). Hybrid MPI-OpenMP paradigm on SMP clusters: MPEG-2 encoder and N-body simulation. arXiv. https://arxiv.org/abs/1211.2292
