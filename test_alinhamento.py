import io
import random
import unittest
import zipfile

import alinhamento as al
from app import app


def rodar(s1, s2, modo, match=2, mismatch=-1, gap=-2):
    return al.executar(s1, s2, match, mismatch, gap, modo)


class TestGlobal(unittest.TestCase):
    def test_exemplo_da_especificacao(self):
        r = rodar("ACGTAC", "ACGTTC", "GLOBAL")
        self.assertEqual(r["score"], 9)
        self.assertEqual(r["seq1_alinhada"], "ACGTAC")
        self.assertEqual(r["seq2_alinhada"], "ACGTTC")
        self.assertEqual(r["linha"], "||||.|")
        self.assertEqual(r["matriz"][0], [0, -2, -4, -6, -8, -10, -12])
        self.assertEqual(r["matriz"][1], [-2, 2, 0, -2, -4, -6, -8])
        self.assertEqual(r["matriz"][6], [-12, -8, -4, 0, 4, 5, 9])
        self.assertEqual(r["traceback"][0], [".", "←", "←", "←", "←", "←", "←"])
        self.assertEqual([r["traceback"][i][0] for i in range(1, 7)], ["↑"] * 6)

    def test_empate_prioriza_diagonal_depois_vertical(self):
        # celula (2,1): diagonal = 1 e vertical = 1 -> deve escolher diagonal
        r = rodar("AA", "A", "GLOBAL", match=2, mismatch=-1, gap=-1)
        self.assertEqual(r["traceback"][2][1], "↖")
        self.assertEqual(r["seq1_alinhada"], "AA")
        self.assertEqual(r["seq2_alinhada"], "-A")

    def test_tamanhos_diferentes(self):
        r = rodar("ACGT", "AC", "GLOBAL")
        self.assertEqual(len(r["seq1_alinhada"]), len(r["seq2_alinhada"]))
        self.assertEqual(r["seq1_alinhada"].replace("-", ""), "ACGT")
        self.assertEqual(r["seq2_alinhada"].replace("-", ""), "AC")


class TestLocal(unittest.TestCase):
    def test_exemplo_classico(self):
        r = rodar("TGTTACGG", "GGTTGACTA", "LOCAL", match=3, mismatch=-3, gap=-2)
        self.assertEqual(r["score"], 13)
        self.assertEqual(r["seq1_alinhada"].replace("-", ""), "GTTAC")
        self.assertEqual(r["seq2_alinhada"].replace("-", ""), "GTTGAC")
        self.assertTrue(all(v >= 0 for linha in r["matriz"] for v in linha))

    def test_sem_similaridade(self):
        r = rodar("AAAA", "CCCC", "LOCAL")
        self.assertEqual(r["score"], 0)
        self.assertEqual(r["seq1_alinhada"], "")
        self.assertEqual(r["tamanho"], 0)

    def test_score_e_o_maior_valor_da_matriz(self):
        r = rodar("ACGTACGT", "TTACGTAA", "LOCAL")
        self.assertEqual(r["score"], max(max(l) for l in r["matriz"]))


class TestValidacao(unittest.TestCase):
    def erros(self, *args):
        with self.assertRaises(al.ErroValidacao) as c:
            al.executar(*args)
        return c.exception.erros

    def test_base_invalida(self):
        e = self.erros("ACGTAXTG", "ACGT", 2, -1, -2, "GLOBAL")
        self.assertIn("base inválida: X", e[0])

    def test_sequencia_vazia(self):
        e = self.erros("", "ACGT", 2, -1, -2, "GLOBAL")
        self.assertIn("vazia", e[0])

    def test_parametro_ausente_e_nao_numerico(self):
        e = self.erros("ACGT", "ACGT", None, "a", -2, "GLOBAL")
        self.assertEqual(len(e), 2)

    def test_parametro_nao_inteiro(self):
        e = self.erros("ACGT", "ACGT", "2.5", -1, -2, "GLOBAL")
        self.assertIn("inteiro", e[0])

    def test_modo_invalido(self):
        self.erros("ACGT", "ACGT", 2, -1, -2, "XYZ")

    def test_minusculas_e_parametros_texto(self):
        r = al.executar("acgt", "acgt", "+2", "-1", "-2", "global")
        self.assertEqual(r["seq1"], "ACGT")
        self.assertEqual(r["score"], 8)

    def test_avisos(self):
        r = al.executar("ACGT", "ACGT", 2, 1, -2, "GLOBAL")
        self.assertEqual(len(r["avisos"]), 1)

    def test_fasta(self):
        self.assertEqual(al.interpretar_texto_sequencias(">Seq1\nACGTAC\n>Seq2\nACGTTC\n"),
                         ("ACGTAC", "ACGTTC"))
        self.assertEqual(al.interpretar_texto_sequencias(">a\nAC\nGT\n>b\nACGT"),
                         ("ACGT", "ACGT"))
        with self.assertRaises(al.ErroValidacao):
            al.interpretar_texto_sequencias(">a\nACGT\n>b\nACGT\n>c\nACGT")
        with self.assertRaises(al.ErroValidacao):
            al.interpretar_texto_sequencias("")
        with self.assertRaises(al.ErroValidacao):
            al.ler_conteudo_arquivo("dados.pdf", b"ACGT")
        with self.assertRaises(al.ErroValidacao):
            al.ler_arquivo("/nao/existe.fasta")


class TestExportacao(unittest.TestCase):
    def test_arquivos(self):
        r = rodar("ACGTAC", "ACGTTC", "GLOBAL")
        arqs = al.gerar_arquivos(r)
        self.assertEqual(arqs["txt"][0], "alinhamento_ACGTAC_ACGTTC.txt")
        self.assertEqual(arqs["pontuacao"][0], "matriz_pontuacao_ACGTAC_ACGTTC.csv")
        txt = arqs["txt"][1].decode()
        self.assertIn("Score final: 9", txt)
        self.assertIn("Seq1: A C G T A C\n      | | | | . |\nSeq2: A C G T T C", txt)
        csv_txt = arqs["alinhamento"][1].decode("utf-8-sig")
        self.assertIn("5,A,.,T", csv_txt)
        z = zipfile.ZipFile(io.BytesIO(al.gerar_zip(r)))
        self.assertEqual(len(z.namelist()), 4)


class TestFlask(unittest.TestCase):
    def setUp(self):
        self.c = app.test_client()
        self.corpo = {"seq1": "ACGTAC", "seq2": "ACGTTC", "modo": "GLOBAL",
                      "match": 2, "mismatch": -1, "gap": -2}

    def test_alinhar(self):
        r = self.c.post("/api/alinhar", json=self.corpo)
        self.assertEqual(r.status_code, 200)
        self.assertEqual(r.get_json()["score"], 9)

    def test_erro_400(self):
        r = self.c.post("/api/alinhar", json={**self.corpo, "seq1": "ACGX"})
        self.assertEqual(r.status_code, 400)
        self.assertTrue(r.get_json()["erros"])

    def test_comparar(self):
        r = self.c.post("/api/comparar", json=self.corpo).get_json()
        self.assertEqual(r["global"]["metodo"], "GLOBAL")
        self.assertEqual(r["local"]["metodo"], "LOCAL")

    def test_upload_fasta(self):
        dados = {"modo": "LOCAL", "match": "2", "mismatch": "-1", "gap": "-2",
                 "arquivo": (io.BytesIO(b">Seq1\nACGTAC\n>Seq2\nACGTTC\n"), "x.fasta")}
        r = self.c.post("/api/alinhar", data=dados, content_type="multipart/form-data")
        self.assertEqual(r.status_code, 200)

    def test_exportar(self):
        for tipo, mime in (("txt", "text/plain"), ("pontuacao", "text/csv"),
                           ("todos", "application/zip")):
            r = self.c.post(f"/api/exportar/{tipo}", json=self.corpo)
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.mimetype, mime)
        self.assertEqual(self.c.post("/api/exportar/xyz", json=self.corpo).status_code, 400)

    def test_pagina_inicial(self):
        r = self.c.get("/")
        self.assertEqual(r.status_code, 200)
        self.assertIn("Passo a passo", r.get_data(as_text=True))


class TestPassos(unittest.TestCase):
    """O passo a passo deve reproduzir exatamente a matriz e o traceback."""

    def _res(self, s1, s2, m, mm, g, modo):
        return al.executar(s1, s2, m, mm, g, modo)

    def test_preenchimento_reproduz_matriz(self):
        for modo in ("GLOBAL", "LOCAL"):
            r = self._res("TGTTACGG", "GGTTGACTA", 3, -2, -2, modo)
            celulas = [p for p in r["passos"]["preenchimento"] if p["tipo"] == "celula"]
            self.assertEqual(len(celulas), 8 * 9)
            for p in celulas:
                self.assertEqual(p["valor"], r["matriz"][p["i"]][p["j"]])
                self.assertEqual(p["dir"], r["traceback"][p["i"]][p["j"]])
                melhor = max(p["diag"], p["cima"], p["esq"])
                self.assertEqual(p["valor"], max(melhor, 0) if modo == "LOCAL" else melhor)

    def test_traceback_reconstroi_alinhamento(self):
        for modo in ("GLOBAL", "LOCAL"):
            r = self._res("TGTTACGG", "GGTTGACTA", 3, -2, -2, modo)
            tb = r["passos"]["traceback"]
            self.assertEqual("".join(x["a"] for x in reversed(tb)), r["seq1_alinhada"])
            self.assertEqual("".join(x["b"] for x in reversed(tb)), r["seq2_alinhada"])
            self.assertEqual(len(tb), r["tamanho"])

    def test_sem_passos_para_matriz_grande(self):
        r = self._res("A" * 150, "C" * 150, 2, -1, -2, "GLOBAL")
        self.assertIsNone(r["passos"])


# ----------------------------------------------------------------------------
# Conferência contra a definição dos algoritmos
# ----------------------------------------------------------------------------
def referencia(s1, s2, match, mismatch, gap, local):
    """
    Implementação escrita direto da definição, separada da do módulo:

      Needleman-Wunsch: F(i,0) = i*gap, F(0,j) = j*gap
                        F(i,j) = max(diag + s(a,b), cima + gap, esq + gap)
      Smith-Waterman:   H(i,0) = H(0,j) = 0
                        H(i,j) = max(0, diag + s(a,b), cima + gap, esq + gap)

    Empate: diagonal > cima > esquerda. No local, a célula que fica <= 0 vira 0
    e não aponta para lugar nenhum. Traceback: global parte de (m,n) e vai até (0,0);
    local parte da primeira célula de maior valor (linha a linha) e para no primeiro 0.
    """
    m, n = len(s1), len(s2)
    H, T = {}, {}
    for i in range(m + 1):
        H[i, 0], T[i, 0] = (0, ".") if local or i == 0 else (i * gap, "↑")
    for j in range(1, n + 1):
        H[0, j], T[0, j] = (0, ".") if local else (j * gap, "←")
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if s1[i - 1] == s2[j - 1] else mismatch
            candidatos = [(H[i - 1, j - 1] + s, "↖"), (H[i - 1, j] + gap, "↑"),
                          (H[i, j - 1] + gap, "←")]
            valor = max(v for v, _ in candidatos)
            seta = next(d for v, d in candidatos if v == valor)  # 1º na ordem de prioridade
            if local and valor <= 0:
                valor, seta = 0, "."
            H[i, j], T[i, j] = valor, seta

    if local:
        maior = max(H.values())
        inicio = next((i, j) for i in range(m + 1) for j in range(n + 1) if H[i, j] == maior)
    else:
        inicio = (m, n)

    i, j = inicio
    a = b = ""
    while (H[i, j] != 0) if local else (i, j) != (0, 0):
        seta = T[i, j]
        a = (s1[i - 1] if seta in "↖↑" else "-") + a
        b = (s2[j - 1] if seta in "↖←" else "-") + b
        i, j = (i - 1 if seta in "↖↑" else i), (j - 1 if seta in "↖←" else j)

    return {
        "matriz": [[H[i, j] for j in range(n + 1)] for i in range(m + 1)],
        "traceback": [[T[i, j] for j in range(n + 1)] for i in range(m + 1)],
        "inicio": list(inicio), "fim": [i, j], "score": H[inicio],
        "seq1_alinhada": a, "seq2_alinhada": b,
    }


class TestDefinicaoDosAlgoritmos(unittest.TestCase):
    CAMPOS = ("matriz", "traceback", "inicio", "fim", "score", "seq1_alinhada", "seq2_alinhada")

    def test_local_calculado_a_mao(self):
        # Seq1 = ACG (linhas), Seq2 = CG (colunas), +2 / -1 / -2
        #        -  C  G
        #   -    0  0  0
        #   A    0  0  0     A×C e A×G: tudo <= 0 -> zera
        #   C    0  2  0     C×C: 0+2 = 2 (↖); C×G: max(-1, -2, 0) = 0 -> zera
        #   G    0  0  4     G×C: max(-1, 0, -2) = 0 -> zera; G×G: 2+2 = 4 (↖)
        r = rodar("ACG", "CG", "LOCAL")
        self.assertEqual(r["matriz"], [[0, 0, 0], [0, 0, 0], [0, 2, 0], [0, 0, 4]])
        self.assertEqual(r["traceback"], [[".", ".", "."], [".", ".", "."],
                                          [".", "↖", "."], [".", ".", "↖"]])
        self.assertEqual((r["inicio"], r["fim"], r["score"]), ([3, 2], [1, 0], 4))
        self.assertEqual((r["seq1_alinhada"], r["seq2_alinhada"]), ("CG", "CG"))
        self.assertEqual((r["regiao_seq1"], r["regiao_seq2"]), ([2, 3], [1, 2]))
        self.assertEqual(r["caminho"], [[2, 1], [3, 2]])

    def test_local_empate_no_maior_valor_comeca_no_primeiro(self):
        # "A" aparece duas vezes em Seq2: as duas células valem 2; vale a 1ª lida linha a linha
        r = rodar("A", "ACA", "LOCAL")
        self.assertEqual(r["matriz"][1], [0, 2, 0, 2])
        self.assertEqual(r["inicio"], [1, 1])

    def test_local_segue_a_definicao_mesmo_com_regras_fora_do_comum(self):
        # Com Gap positivo o Smith-Waterman continua seguindo a fórmula:
        # H(1,1) = max(0, 0+1, 0+1, 0+1) = 1 e o empate fica com a diagonal.
        r = rodar("C", "A", "LOCAL", match=2, mismatch=1, gap=1)
        self.assertEqual(r["matriz"], [[0, 0], [0, 1]])
        self.assertEqual(r["traceback"][1][1], "↖")
        self.assertEqual(r["score"], 1)

    def _compara_com_referencia(self, modo, casos=1500, semente=2024):
        rnd = random.Random(semente)
        regras = [(m, x, g) for m in (-1, 0, 1, 2, 3, 5) for x in (-3, -1, 0, 1, 2)
                  for g in (-4, -2, -1, 0, 1, 2)]
        for _ in range(casos):
            s1 = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(1, 15)))
            s2 = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(1, 15)))
            m, x, g = rnd.choice(regras)
            obtido = rodar(s1, s2, modo, match=m, mismatch=x, gap=g)
            esperado = referencia(s1, s2, m, x, g, modo == "LOCAL")
            for campo in self.CAMPOS:
                self.assertEqual(obtido[campo], esperado[campo],
                                 f"{modo} {s1} x {s2} ({m}/{x}/{g}): campo {campo}")
            self._confere_passos(obtido, modo == "LOCAL")

    def _confere_passos(self, r, local):
        """Os dados do passo a passo (usados na tela) mostram o que o algoritmo fez."""
        H, T, s1, s2 = r["matriz"], r["traceback"], r["seq1"], r["seq2"]
        celulas = [p for p in r["passos"]["preenchimento"] if p["tipo"] == "celula"]
        self.assertEqual([(p["i"], p["j"]) for p in celulas],
                         [(i, j) for i in range(1, len(s1) + 1) for j in range(1, len(s2) + 1)])
        for p in celulas:
            i, j = p["i"], p["j"]
            s = r["match"] if s1[i - 1] == s2[j - 1] else r["mismatch"]
            self.assertEqual((p["diag"], p["cima"], p["esq"]),
                             (H[i - 1][j - 1] + s, H[i - 1][j] + r["gap"], H[i][j - 1] + r["gap"]))
            self.assertEqual((p["valor"], p["dir"]), (H[i][j], T[i][j]))
            self.assertEqual(p["zerou"], local and max(p["diag"], p["cima"], p["esq"]) <= 0)
        tb = r["passos"]["traceback"]
        self.assertEqual("".join(p["a"] for p in reversed(tb)), r["seq1_alinhada"])
        self.assertEqual("".join(p["b"] for p in reversed(tb)), r["seq2_alinhada"])
        self.assertEqual(r["passos"]["parada"]["valor"], H[r["fim"][0]][r["fim"][1]])
        if local:
            self.assertEqual(r["passos"]["parada"]["valor"], 0)

    def test_local_igual_a_definicao_em_casos_aleatorios(self):
        self._compara_com_referencia("LOCAL")

    def test_global_igual_a_definicao_em_casos_aleatorios(self):
        self._compara_com_referencia("GLOBAL")


if __name__ == "__main__":
    unittest.main()
