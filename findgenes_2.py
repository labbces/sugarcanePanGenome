import os
import time
import argparse
import logging
from multiprocessing import Pool, cpu_count

def normalize_chromosome(chrom):
    """Padroniza os nomes dos cromossomos para uma comparação mais robusta."""
    chrom = chrom.replace("Chr", "ChrA")  # Ajusta nomes que podem estar abreviados
    chrom = chrom.replace("ChrA0", "ChrA")  # Remove zeros à esquerda
    return chrom

def parse_distritype(file_path):
    """Parse the distritype.txt file into a list of dictionaries."""
    with open(file_path, 'r') as file:
        lines = file.read().strip().split('\n\n')
    parsed_data = []
    for block in lines:
        entry = {}
        for line in block.split('\n'):
            if ': ' in line:
                key, value = line.split(': ', 1)
                entry[key] = value
            else:
                logging.warning(f"Skipping line with unexpected format: {line}")
        if "Chromosome" in entry:
            entry["Chromosome"] = normalize_chromosome(entry["Chromosome"])
        parsed_data.append(entry)
    return parsed_data

def filter_by_haplotype_and_chromosome(data, haplotype):
    """Filter the parsed data to include only entries with the specified haplotype."""
    haplotype = haplotype.strip().lower()
    return [entry for entry in data if entry['Haplotype'].lower() == haplotype]

def process_line(line, data):
    """Process a single line of the GFF3 file and extract genes that overlap with the specified ranges."""
    if line.startswith('#'):
        return None

    columns = line.strip().split('\t')
    if len(columns) < 9:
        logging.warning(f"Skipping line with unexpected format: {line.strip()}")
        return None

    try:
        seqid = normalize_chromosome(columns[0])  # Normaliza cromossomo
        feature_type = columns[2]
        start, end = int(columns[3]), int(columns[4])
    except ValueError as ve:
        logging.error(f"Error parsing columns in line: {line.strip()} - {ve}")
        return None

    # Verifica se é um gene
    if feature_type != 'gene':
        return None

    for entry in data:
        try:
            line_start, line_end = map(int, entry['Line'].split())

            # Correção: Incluir genes que se sobrepõem ao intervalo
            if seqid == entry['Chromosome'] and not (end < line_start or start > line_end):
                category = entry['Max Size Distribution'].split(', ')[0]
                return {
                    'file': line,
                    'category': category,
                    'chromosome': seqid
                }
        except ValueError as ve:
            logging.error(f"Error parsing Line field in distritype entry: {entry['Line']} - {ve}")
            continue

    return None

def extract_genes(gff_file, data, num_cpus):
    """Extract genes from the GFF3 file that overlap with specified ranges."""
    with open(gff_file, 'r') as file:
        lines = file.readlines()

    with Pool(num_cpus) as pool:
        results = pool.starmap(process_line, [(line, data) for line in lines])

    return [result for result in results if result is not None]

def write_results(genes, output_dir, haplotype):
    """Write the extracted genes to separate files based on their categories and chromosomes."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    category_files = {}
    gene_counts = {}

    for gene in genes:
        category = gene['category']
        chromosome = gene['chromosome']
        key = f"{category}_{chromosome}"

        if key not in category_files:
            category_files[key] = open(os.path.join(output_dir, f'{haplotype}_genes_{category.replace(" ", "_")}_{chromosome}.txt'), 'w')
            gene_counts[key] = 0

        category_files[key].write(gene['file'])
        gene_counts[key] += 1

    # Fechar os arquivos
    for file in category_files.values():
        file.close()

    # Exibir e registrar o número de genes por categoria
    for category, count in gene_counts.items():
        print(f'Arquivo gerado para {category}: {count} genes encontrados')
        logging.info(f'Arquivo gerado para {category}: {count} genes encontrados')

def main():
    """Main function to parse command-line arguments and execute the workflow."""
    parser = argparse.ArgumentParser(description='Process GFF3 files based on distritype.txt')
    parser.add_argument('--haplotype', required=True, help='Haplotype to filter')
    parser.add_argument('--distritype', required=True, help='Path to the distritype.txt file')
    parser.add_argument('--gff', required=True, help='Path to the GFF3 file')
    parser.add_argument('--output', required=True, help='Directory to save output files')
    parser.add_argument('--cpus', type=int, default=cpu_count(), help='Number of CPUs to use for parallel processing')
    args = parser.parse_args()

    logging.basicConfig(filename='process2.log', level=logging.INFO)
    start_time = time.time()

    try:
        # Parse distritype.txt
        distritype_data = parse_distritype(args.distritype)
        logging.info(f'Total entries processed from distritype.txt: {len(distritype_data)}')

        # Filter data by haplotype
        filtered_data = filter_by_haplotype_and_chromosome(distritype_data, args.haplotype)

        # Extract genes
        genes = extract_genes(args.gff, filtered_data, args.cpus)

        # Write results
        write_results(genes, args.output, args.haplotype)

        end_time = time.time()
        logging.info(f'Processing completed in {end_time - start_time:.2f} seconds')
        logging.info(f'Total genes found: {len(genes)}')

        # Exibir total de genes encontrados
        print(f'Total de genes extraídos: {len(genes)}')

    except Exception as e:
        logging.error(f'Error: {e}')
        print(f'Erro: {e}')

if __name__ == '__main__':
    main()

