"""
Alinhamento de duas sequências de DNA por programação dinâmica.

  - GLOBAL: Needleman-Wunsch
  - LOCAL : Smith-Waterman

Este módulo contém apenas a lógica (validação, algoritmos, traceback,
formatação e geração dos arquivos de exportação). Não depende do Flask.
"""
import csv
import io
import math
import os
import re
import zipfile
from dataclasses import dataclass, field

# ----------------------------------------------------------------------------
# Constantes
# ----------------------------------------------------------------------------
BASES_VALIDAS = set("ACGT")
GAP_CHAR = "-"

# Direções guardadas na matriz de traceback
DIAGONAL = "↖"      # base de Seq1 alinhada com base de Seq2
VERTICAL = "↑"      # base de Seq1 alinhada com gap em Seq2
HORIZONTAL = "←"    # gap em Seq1 alinhado com base de Seq2
INICIO = "."        # origem (global) ou célula zerada (local)

NOMES_METODO = {
    "GLOBAL": "Global (Needleman-Wunsch)",
    "LOCAL": "Local (Smith-Waterman)",
}
EXTENSOES_ACEITAS = (".txt", ".fasta", ".fa")


# ----------------------------------------------------------------------------
# Validação
# ----------------------------------------------------------------------------
class ErroValidacao(Exception):
    """Erro de entrada. `erros` guarda todas as mensagens encontradas."""

    def __init__(self, erros):
        if isinstance(erros, str):
            erros = [erros]
        self.erros = list(erros)
        super().__init__("; ".join(self.erros))


@dataclass
class Entrada:
    seq1: str
    seq2: str
    match: int
    mismatch: int
    gap: int
    modo: str
    avisos: list = field(default_factory=list)


def normalizar_sequencia(seq, nome):
    """Retorna (sequência_maiúscula, lista_de_erros)."""
    if not isinstance(seq, str) or seq.strip() == "":
        return None, [f"{nome}: sequência vazia."]
    seq = seq.strip().upper()
    erros, ja_reportados = [], set()
    for pos, c in enumerate(seq, start=1):
        if c not in BASES_VALIDAS and c not in ja_reportados:
            ja_reportados.add(c)
            erros.append(f"ERRO – base inválida: {c} ({nome}, posição {pos})")
    return (None if erros else seq), erros


def validar_parametro(valor, nome):
    """Retorna (inteiro, mensagem_de_erro). Aceita int ou texto como '+2', '-1'."""
    if valor is None or (isinstance(valor, str) and valor.strip() == ""):
        return None, f"Parâmetro ausente: {nome}."
    if isinstance(valor, bool):
        return None, f"Parâmetro não numérico: {nome}."
    if isinstance(valor, int):
        return valor, None
    if isinstance(valor, float):
        if math.isfinite(valor) and valor.is_integer():
            return int(valor), None
        return None, f"Parâmetro {nome} deve ser um número inteiro."
    if isinstance(valor, str):
        s = valor.strip()
        if re.fullmatch(r"[+-]?\d+", s):
            return int(s), None
        try:
            f = float(s)
        except ValueError:
            return None, f"Parâmetro não numérico: {nome} = '{valor}'."
        if not math.isfinite(f):
            return None, f"Parâmetro não numérico: {nome} = '{valor}'."
        return None, f"Parâmetro {nome} deve ser um número inteiro."
    return None, f"Parâmetro não numérico: {nome}."


def validar_entrada(seq1, seq2, match, mismatch, gap, modo):
    """Valida tudo antes de processar. Levanta ErroValidacao com todos os erros."""
    erros = []

    s1, e1 = normalizar_sequencia(seq1, "Seq1")
    s2, e2 = normalizar_sequencia(seq2, "Seq2")
    erros += e1 + e2

    valores = {}
    for nome, bruto in (("Match", match), ("Mismatch", mismatch), ("Gap", gap)):
        v, e = validar_parametro(bruto, nome)
        valores[nome] = v
        if e:
            erros.append(e)

    modo_norm = modo.strip().upper() if isinstance(modo, str) else None
    if modo_norm not in NOMES_METODO:
        erros.append("Tipo de alinhamento inválido ou ausente: use LOCAL ou GLOBAL.")

    if erros:
        raise ErroValidacao(erros)

    avisos = []
    if valores["Match"] <= 0:
        avisos.append("Recomendado: Match > 0.")
    if valores["Mismatch"] >= 0:
        avisos.append("Recomendado: Mismatch < 0.")
    if valores["Gap"] >= 0:
        avisos.append("Recomendado: Gap < 0.")

    return Entrada(s1, s2, valores["Match"], valores["Mismatch"], valores["Gap"],
                   modo_norm, avisos)


# ----------------------------------------------------------------------------
# Leitura de arquivos (.txt / .fasta)
# ----------------------------------------------------------------------------
def interpretar_texto_sequencias(texto):
    """
    Extrai exatamente duas sequências de um texto.
    Aceita FASTA (>Seq1 / >Seq2, sequência podendo ocupar várias linhas)
    ou, sem cabeçalhos '>', uma sequência por linha.
    """
    linhas = [l.strip() for l in texto.splitlines() if l.strip()]
    if not linhas:
        raise ErroValidacao("Arquivo em formato incompatível: nenhuma sequência encontrada.")

    if any(l.startswith(">") for l in linhas):
        if not linhas[0].startswith(">"):
            raise ErroValidacao(
                "Arquivo em formato incompatível: há conteúdo antes do primeiro cabeçalho '>'.")
        registros = []
        for l in linhas:
            if l.startswith(">"):
                registros.append([])
            else:
                registros[-1].append(l)
        seqs = ["".join(r) for r in registros]
    else:
        seqs = linhas

    if len(seqs) != 2:
        raise ErroValidacao(
            f"Quantidade de sequências diferente de duas: encontradas {len(seqs)}.")
    return seqs[0], seqs[1]


def ler_conteudo_arquivo(nome_arquivo, dados: bytes):
    """Usado para arquivos enviados (upload). Retorna (seq1, seq2) ainda sem validar bases."""
    if os.path.splitext(nome_arquivo or "")[1].lower() not in EXTENSOES_ACEITAS:
        raise ErroValidacao("Arquivo em formato incompatível: use .txt ou .fasta.")
    try:
        texto = dados.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise ErroValidacao("Arquivo ilegível: não foi possível decodificar como texto.")
    return interpretar_texto_sequencias(texto)


def ler_arquivo(caminho):
    """Lê um arquivo no disco. Retorna (seq1, seq2)."""
    if not os.path.isfile(caminho):
        raise ErroValidacao(f"Arquivo inexistente: {caminho}")
    try:
        with open(caminho, "rb") as f:
            dados = f.read()
    except OSError:
        raise ErroValidacao(f"Arquivo ilegível: {caminho}")
    return ler_conteudo_arquivo(caminho, dados)


# ----------------------------------------------------------------------------
# Algoritmos
# ----------------------------------------------------------------------------
def _pontuar(a, b, match, mismatch):
    return match if a == b else mismatch


def preencher_matriz(seq1, seq2, match, mismatch, gap, local):
    """
    Constrói a matriz de pontuação H e a matriz de direções T.
    Linhas = Seq1 (i), colunas = Seq2 (j). Ambas têm (m+1) x (n+1).

    Prioridade em empates: diagonal > vertical > horizontal.
    """
    m, n = len(seq1), len(seq2)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    T = [[INICIO] * (n + 1) for _ in range(m + 1)]

    if not local:  # Needleman-Wunsch: gaps acumulados na borda
        for i in range(1, m + 1):
            H[i][0] = i * gap
            T[i][0] = VERTICAL
        for j in range(1, n + 1):
            H[0][j] = j * gap
            T[0][j] = HORIZONTAL

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i - 1][j - 1] + _pontuar(seq1[i - 1], seq2[j - 1], match, mismatch)
            cima = H[i - 1][j] + gap
            esq = H[i][j - 1] + gap
            melhor = max(diag, cima, esq)

            if local and melhor <= 0:      # max(0, ...) do Smith-Waterman
                H[i][j] = 0
                T[i][j] = INICIO
                continue

            H[i][j] = melhor
            if diag == melhor:
                T[i][j] = DIAGONAL
            elif cima == melhor:
                T[i][j] = VERTICAL
            else:
                T[i][j] = HORIZONTAL
    return H, T


def _ponto_inicial(H, local):
    """Global: canto inferior direito. Local: maior valor (empate: 1º em ordem linha a linha)."""
    m, n = len(H) - 1, len(H[0]) - 1
    if not local:
        return m, n
    melhor, pos = H[0][0], (0, 0)
    for i in range(m + 1):
        for j in range(n + 1):
            if H[i][j] > melhor:
                melhor, pos = H[i][j], (i, j)
    return pos


def _traceback(seq1, seq2, H, T, inicio, local):
    """
    Retorna (seq1_alinhada, seq2_alinhada, caminho, fim).
    `caminho` lista as células que geraram uma coluna do alinhamento (para destaque
    na matriz); `fim` é a célula onde o traceback terminou.
    """
    i, j = inicio
    a, b, caminho = [], [], []
    while True:
        if local and H[i][j] == 0:
            break
        if not local and i == 0 and j == 0:
            break
        caminho.append((i, j))
        d = T[i][j]
        if d == DIAGONAL:
            a.append(seq1[i - 1]); b.append(seq2[j - 1]); i -= 1; j -= 1
        elif d == VERTICAL:
            a.append(seq1[i - 1]); b.append(GAP_CHAR); i -= 1
        else:  # HORIZONTAL
            a.append(GAP_CHAR); b.append(seq2[j - 1]); j -= 1
    caminho.reverse()
    return "".join(reversed(a)), "".join(reversed(b)), caminho, (i, j)


def _linha_marcadores(a, b):
    out = []
    for x, y in zip(a, b):
        if x == GAP_CHAR or y == GAP_CHAR:
            out.append(" ")
        elif x == y:
            out.append("|")
        else:
            out.append(".")
    return "".join(out)


# Limite de células (m x n) para gerar o passo a passo detalhado na resposta.
LIMITE_PASSOS = 10000


def gerar_passos(seq1, seq2, H, T, caminho, inicio, fim, match, mismatch, gap, local):
    """
    Reconstrói, a partir das matrizes prontas, cada decisão do algoritmo:
      - "preenchimento": uma entrada por célula, na ordem em que é calculada
        (bordas incluídas), com os três candidatos e a direção escolhida;
      - "traceback": uma entrada por movimento, do ponto inicial até a parada.
    """
    m, n = len(seq1), len(seq2)
    preench = []

    if local:
        preench.append({"tipo": "init", "i": 0, "j": 0, "valor": 0, "bordas": True})
    else:
        preench.append({"tipo": "init", "i": 0, "j": 0, "valor": 0, "bordas": False})
        for j in range(1, n + 1):
            preench.append({"tipo": "init", "i": 0, "j": j, "valor": H[0][j],
                            "dir": HORIZONTAL})
        for i in range(1, m + 1):
            preench.append({"tipo": "init", "i": i, "j": 0, "valor": H[i][0],
                            "dir": VERTICAL})

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            a, b = seq1[i - 1], seq2[j - 1]
            sc = _pontuar(a, b, match, mismatch)
            diag = H[i - 1][j - 1] + sc
            cima = H[i - 1][j] + gap
            esq = H[i][j - 1] + gap
            preench.append({
                "tipo": "celula", "i": i, "j": j, "a": a, "b": b,
                "score": sc, "igual": a == b,
                "diag": diag, "cima": cima, "esq": esq,
                "diag_de": H[i - 1][j - 1], "cima_de": H[i - 1][j],
                "esq_de": H[i][j - 1],
                "valor": H[i][j], "dir": T[i][j],
                # só no local: o max(0, ...) venceu os três candidatos
                "zerou": bool(local and max(diag, cima, esq) <= 0),
            })

    tb = []
    for (i, j) in reversed(caminho):
        d = T[i][j]
        tb.append({
            "i": i, "j": j, "valor": H[i][j], "dir": d,
            "a": seq1[i - 1] if d in (DIAGONAL, VERTICAL) else GAP_CHAR,
            "b": seq2[j - 1] if d in (DIAGONAL, HORIZONTAL) else GAP_CHAR,
        })
    return {"preenchimento": preench, "traceback": tb,
            "parada": {"i": fim[0], "j": fim[1], "valor": H[fim[0]][fim[1]]}}



def alinhar(entrada: Entrada):
    """Executa o algoritmo escolhido em uma Entrada já validada."""
    local = entrada.modo == "LOCAL"
    s1, s2 = entrada.seq1, entrada.seq2

    H, T = preencher_matriz(s1, s2, entrada.match, entrada.mismatch, entrada.gap, local)
    inicio = _ponto_inicial(H, local)
    a, b, caminho, fim = _traceback(s1, s2, H, T, inicio, local)
    linha = _linha_marcadores(a, b)

    passos = None
    if len(s1) * len(s2) <= LIMITE_PASSOS:
        passos = gerar_passos(s1, s2, H, T, caminho, inicio, fim,
                              entrada.match, entrada.mismatch, entrada.gap, local)

    return {
        "metodo": entrada.modo,
        "nome_metodo": NOMES_METODO[entrada.modo],
        "seq1": s1,
        "seq2": s2,
        "match": entrada.match,
        "mismatch": entrada.mismatch,
        "gap": entrada.gap,
        "matriz": H,
        "traceback": T,
        "inicio": list(inicio),          # célula onde o traceback começou
        "fim": list(fim),                # célula onde terminou
        "caminho": [list(c) for c in caminho],
        "seq1_alinhada": a,
        "linha": linha,
        "seq2_alinhada": b,
        "score": H[inicio[0]][inicio[1]],
        "tamanho": len(a),
        "matches": linha.count("|"),
        "mismatches": linha.count("."),
        "gaps": a.count(GAP_CHAR) + b.count(GAP_CHAR),
        # trechos (1-based, inclusivos) de cada sequência usados no alinhamento
        "regiao_seq1": [fim[0] + 1, inicio[0]] if inicio[0] > fim[0] else None,
        "regiao_seq2": [fim[1] + 1, inicio[1]] if inicio[1] > fim[1] else None,
        "avisos": entrada.avisos,
        "passos": passos,               # None quando m*n > LIMITE_PASSOS
    }


def executar(seq1, seq2, match, mismatch, gap, modo):
    """Valida a entrada e só então processa. Ponto de entrada principal."""
    return alinhar(validar_entrada(seq1, seq2, match, mismatch, gap, modo))


# ----------------------------------------------------------------------------
# Formatação de texto
# ----------------------------------------------------------------------------
def formatar_alinhamento(res):
    """Três linhas: Seq1 / marcadores / Seq2, no estilo da especificação."""
    a = " ".join(res["seq1_alinhada"])
    m = " ".join(res["linha"])
    b = " ".join(res["seq2_alinhada"])
    return f"Seq1: {a}\n      {m}\nSeq2: {b}"


def formatar_matriz_texto(matriz, seq1, seq2):
    """Matriz (ou traceback) em texto alinhado, útil para modo console."""
    linhas = [["", "-", *seq2]]
    for i, row in enumerate(matriz):
        linhas.append([("-" if i == 0 else seq1[i - 1]), *[str(v) for v in row]])
    larg = max(len(c) for l in linhas for c in l)
    return "\n".join(" ".join(c.rjust(larg) for c in l) for l in linhas)


def gerar_txt(res):
    """Conteúdo do arquivo de texto (mesmo formato da tela de exportação)."""
    return (
        f"Algoritmo: {res['nome_metodo']}\n"
        f"Sequência 1: {res['seq1']}\n"
        f"Sequência 2: {res['seq2']}\n"
        f"Match:     {res['match']}\n"
        f"Mismatch: {res['mismatch']}\n"
        f"Gap: {res['gap']}\n"
        f"Score final: {res['score']}\n"
        f"Tamanho do alinhamento: {res['tamanho']}\n"
        "\n"
        "Alinhamento final:\n"
        f"{formatar_alinhamento(res)}\n"
        "\n"
        "Matriz de pontuação:\n(ver arquivo CSV correspondente)\n"
        "\n"
        "Matriz de traceback:\n(ver arquivo CSV correspondente)\n"
    )


# ----------------------------------------------------------------------------
# CSV / arquivos de exportação
# ----------------------------------------------------------------------------
def _csv(linhas):
    buf = io.StringIO()
    csv.writer(buf, lineterminator="\n").writerows(linhas)
    return buf.getvalue()


def matriz_para_csv(matriz, seq1, seq2):
    """Serve tanto para a matriz de pontuação quanto para a de traceback."""
    linhas = [["", "-", *seq2]]
    for i, row in enumerate(matriz):
        linhas.append([("-" if i == 0 else seq1[i - 1]), *row])
    return _csv(linhas)


def alinhamento_para_csv(res):
    linhas = [["Posição", "Seq1", "Linha", "Seq2"]]
    for k, (x, l, y) in enumerate(
            zip(res["seq1_alinhada"], res["linha"], res["seq2_alinhada"]), start=1):
        linhas.append([k, x, l.strip(), y])
    return _csv(linhas)


def nomes_arquivos(res):
    def cortar(s):
        return s[:30]
    sufixo = f"{cortar(res['seq1'])}_{cortar(res['seq2'])}"
    return {
        "txt": f"alinhamento_{sufixo}.txt",
        "pontuacao": f"matriz_pontuacao_{sufixo}.csv",
        "traceback": f"matriz_traceback_{sufixo}.csv",
        "alinhamento": f"alinhamento_{sufixo}.csv",
    }


def gerar_arquivos(res):
    """Retorna {tipo: (nome_do_arquivo, bytes)} para os quatro arquivos mínimos."""
    nomes = nomes_arquivos(res)
    conteudos = {
        "txt": gerar_txt(res),
        "pontuacao": matriz_para_csv(res["matriz"], res["seq1"], res["seq2"]),
        "traceback": matriz_para_csv(res["traceback"], res["seq1"], res["seq2"]),
        "alinhamento": alinhamento_para_csv(res),
    }
    # utf-8-sig nos CSV para o Excel exibir corretamente as setas
    return {
        tipo: (nomes[tipo],
               texto.encode("utf-8-sig" if tipo != "txt" else "utf-8"))
        for tipo, texto in conteudos.items()
    }


def gerar_zip(res):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for nome, dados in gerar_arquivos(res).values():
            z.writestr(nome, dados)
    return buf.getvalue()
