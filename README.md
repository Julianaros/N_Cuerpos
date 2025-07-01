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
