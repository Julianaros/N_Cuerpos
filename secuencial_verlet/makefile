# Makefile para limpiar archivos generados

# Archivos a eliminar
PNGS = $(wildcard *.png)
CSVS = $(wildcard *.csv)

# Elimina solo archivos PNG
clean-png:
	@echo "Eliminando archivos PNG..."
	@rm -f $(PNGS)

# Elimina solo archivos CSV
clean-csv:
	@echo "Eliminando archivos CSV..."
	@rm -f $(CSVS)

# Elimina todos los archivos .png y .csv
clean-all: clean-png clean-csv
	@echo "Limpieza completa realizada."

.PHONY: clean-png clean-csv clean-all
