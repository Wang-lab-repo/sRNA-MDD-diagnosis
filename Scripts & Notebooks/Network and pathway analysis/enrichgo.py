import argparse
import scipy.stats as stats

# 使用 argparse 获取命令行输入
def parse_args():
    parser = argparse.ArgumentParser(description="Process gene ontology data and compute Fisher's exact test.")
    parser.add_argument('goa_file', type=str, help="Path to the goa_human.gaf file")
    parser.add_argument('obo_file', type=str, help="Path to the go.obo file")
    parser.add_argument('fasta_file', type=str, help="Path to the uniprot-proteome fasta file")
    parser.add_argument('gene_cor_file', type=str, help="Path to the gene-cor1.txt file")
    parser.add_argument('output_file', type=str, help="Path to the output file (go-cor.txt)")
    return parser.parse_args()

def main():
    # 获取命令行参数
    args = parse_args()

    # 数据结构初始化
    go = {}
    gog = {}
    
    # 读取 goa_human.gaf 文件
    with open(args.goa_file, 'r') as f:
        for l in f:
            sp = l.rstrip().split('\t')
            if 'UniProtKB' in sp:
                if sp[4] in go:
                    go[sp[4]] = go[sp[4]] + [sp[1]]
                else:
                    go[sp[4]] = [sp[1]]
                gog[sp[1]] = ''

    goid = {}
    gotp = {}
    gid = ''
    name = ''
    tp = ''
    
    # 读取 go.obo 文件
    with open(args.obo_file, 'r') as f:
        for l in f:
            l = l.rstrip()
            if l == '[Term]':
                if gid != '' and name != '' and tp != '':
                    goid[gid] = name
                    gotp[gid] = tp
                gid = ''
                name = ''
                tp = ''
            if l.startswith('id: GO:'):
                gid = l.split('id: ')[1]
            if l.startswith('name: '):
                name = l.split('name: ')[1][0].upper() + l.split('name: ')[1][1:]
            if l.startswith('namespace: '):
                tp = l.split('namespace: ')[1]

    NN = {}
    bgene = {}
    geneid = {}
    
    # 读取 uniprot-proteome fasta 文件
    with open(args.fasta_file, 'r') as f:
        for l in f:
            l = l.rstrip()
            if l.startswith('>'):
                id = l.split('|')[1]
                if 'GN=' in l:
                    gene = l.split('GN=')[1].split(' ')[0]
                    geneid[gene] = id
                if id in gog:
                    NN[id] = ''
                bgene[id] = ''

    uid = {}
    MM = {}
    gene = {}
    
    # 读取 gene-cor1.txt 文件
    with open(args.gene_cor_file, 'r', encoding='utf-8') as f:
        for l in f:
            sp = l.rstrip().split('\t')
            sp[1] = sp[1].upper()
            if sp[1] not in geneid:
                continue
            sp[0] = geneid[sp[1]]
            if sp[0] not in bgene:
                bgene[sp[0]] = ''
                if sp[0] in gog:
                    NN[sp[0]] = ''
            if sp[0] in gog:
                MM[sp[0]] = ''

            gene[sp[0]] = ''
            uid[sp[0]] = sp[1]

    print(len(gene), len(bgene), len(goid), len(go), len(MM), len(NN))

    # 写入输出文件
    with open(args.output_file, 'w') as w:
        for gg in go:
            mm, nn = {}, {}
            mm1 = {}
            for gx in MM:
                if gx in go[gg]:
                    mm[gx] = ''
                    mm1[uid[gx]] = ''
            for gx in NN:
                if gx in go[gg]:
                    nn[gx] = ''
            m = len(mm)
            M = len(MM)
            n = len(nn)
            N = len(NN)

            p = stats.fisher_exact([[m, M - m], [n - m, N - n - M + m]])
            if M == 0 or N == 0:
                continue

            w.write(f"{gg}\t{goid[gg]}\t{gotp[gg]}\t{m}\t{M}\t{m / M}\t{n}\t{N}\t{n / N}\t{(m / M) / (n / N)}\t{p[1]}\t{', '.join(mm1)}\n")
            print(p)

if __name__ == "__main__":
    main()
