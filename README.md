# complicator
Calculadora de indices de complejidad para tratamientos VMAT e IMRT

## Interfaz en línea de comandos
Para procesar en línea de comandos un fichero RTPlan y obtener los índices de complejidad
hay que hacer

```
python cli.py fichero_entrada`
```

Si además queremos que el programa genere un fichero debug.txt con información para verificar 
el cálculo, hay que hacer

```
python cli.py --debug fichero_entrada
```

La lista completa de opciones se puede obtener con

```
python cli.py --help
```
