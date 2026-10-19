"""
Aplicação Flask do alinhamento de duas sequências: interface web (GET /) + API.

Endpoints (todos POST):
  /api/alinhar            -> executa o modo escolhido (LOCAL ou GLOBAL)
  /api/comparar           -> executa LOCAL e GLOBAL para o mesmo par de sequências
  /api/exportar/<tipo>    -> baixa arquivo: txt | pontuacao | traceback | alinhamento | todos (zip)

Corpo da requisição:
  JSON:       {"seq1": "ACGTAC", "seq2": "ACGTTC", "modo": "GLOBAL",
               "match": 2, "mismatch": -1, "gap": -2}
  ou multipart/form-data com os mesmos campos; o campo "arquivo" (.txt/.fasta)
  substitui seq1/seq2 quando enviado.

Erros de validação retornam HTTP 400: {"erros": ["...", "..."]}
"""
import io

from flask import Flask, jsonify, render_template, request, send_file

import alinhamento as al

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024  # 5 MB
app.json.ensure_ascii = False                       # mantém acentos e setas no JSON


def _ler_requisicao():
    if request.is_json:
        dados = request.get_json(silent=True)
        if not isinstance(dados, dict):
            raise al.ErroValidacao("Corpo JSON inválido.")
        arquivo = None
    else:
        dados = request.form.to_dict()
        arquivo = request.files.get("arquivo")

    seq1, seq2 = dados.get("seq1"), dados.get("seq2")
    if arquivo is not None and arquivo.filename:
        seq1, seq2 = al.ler_conteudo_arquivo(arquivo.filename, arquivo.read())

    return {
        "seq1": seq1, "seq2": seq2,
        "match": dados.get("match"),
        "mismatch": dados.get("mismatch"),
        "gap": dados.get("gap"),
        "modo": dados.get("modo"),
    }


@app.errorhandler(al.ErroValidacao)
def _erro_validacao(e):
    return jsonify({"erros": e.erros}), 400


@app.errorhandler(413)
def _arquivo_grande(_):
    return jsonify({"erros": ["Arquivo muito grande (limite de 5 MB)."]}), 413


@app.get("/")
def index():
    return render_template("index.html")


@app.post("/api/alinhar")
def alinhar():
    return jsonify(al.executar(**_ler_requisicao()))


@app.post("/api/comparar")
def comparar():
    args = _ler_requisicao()
    args.pop("modo")
    return jsonify({
        "global": al.executar(**args, modo="GLOBAL"),
        "local": al.executar(**args, modo="LOCAL"),
    })


@app.post("/api/exportar/<tipo>")
def exportar(tipo):
    if tipo not in ("txt", "pontuacao", "traceback", "alinhamento", "todos"):
        raise al.ErroValidacao("Tipo de exportação inválido.")
    res = al.executar(**_ler_requisicao())

    if tipo == "todos":
        nome = al.nomes_arquivos(res)["txt"].replace("alinhamento_", "alinhamento_completo_")
        return send_file(io.BytesIO(al.gerar_zip(res)), mimetype="application/zip",
                         as_attachment=True, download_name=nome.rsplit(".", 1)[0] + ".zip")

    nome, dados = al.gerar_arquivos(res)[tipo]
    mime = "text/plain" if tipo == "txt" else "text/csv"
    return send_file(io.BytesIO(dados), mimetype=mime, as_attachment=True, download_name=nome)


if __name__ == "__main__":
    app.run(debug=True)
