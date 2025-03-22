% Huffman 
Alice = fopen('alice_in_wonderland.txt', 'r'); % Abre o arquivo para leitura
Dados = fread(Alice, '*char'); % Lê o conteúdo do arquivo como caracteres
fclose(Alice); % Fecha o arquivo
Total = length(Dados); % Calcula o tamanho total dos dados

% Construir a tabela de frequência
i = 1;
B = [];
while i <= Total
    if JaExiste(Dados(i), B)
        B(end+1, 1) = double(Dados(i));  % Adiciona o caracter ao array B
        B(end, 2) = double(sum(Dados == Dados(i))); % Conta a frequência do caracter
    end
    i = i + 1;
end

ordCar = sortrows(B, -2); % Ordena os caracteres pela frequência em ordem decrescente
Simbolos = char(ordCar(:, 1)); % Obtém os símbolos
ordCar(:, 2) = ordCar(:, 2) / Total; % Normaliza as frequências
numSimb = numel(ordCar(:, 1)); % Número de símbolos
ordSimb = transpose(ordCar(:, 1)); % Transpõe os símbolos
Nos = num2cell(char(ordSimb)); % Converte os símbolos em células
ordProb = transpose(ordCar(:, 2)); % Transpõe as probabilidades

for j = 1:numSimb-1
    [~, indice] = min(ordProb); % Encontra o índice do menor valor de probabilidade
    prob1 = ordProb(indice); % Probabilidade do menor valor
    no1 = Nos{indice}; % Nó do menor valor
    ordProb(indice) = []; % Remove o menor valor de probabilidade
    Nos(indice) = []; % Remove o menor nó

    [~, indice] = min(ordProb); % Repete para o segundo menor valor
    prob2 = ordProb(indice); % Probabilidade do segundo menor valor
    no2 = Nos{indice}; % Nó do segundo menor valor
    ordProb(indice) = []; % Remove o segundo menor valor de probabilidade
    Nos(indice) = []; % Remove o segundo menor nó

    novoNo = {no1, no2}; % Cria um novo nó combinando os dois menores nós
    novaProb = prob1 + prob2; % Soma as probabilidades

    id = find(ordProb >= novaProb, 1, 'last'); % Encontra a posição onde inserir o novo nó

    if isempty(id)
        id = numel(ordProb) + 1; % Se não encontrar, insere no final
    end

    ordProb = [ordProb(1:id-1), novaProb, ordProb(id:end)]; % Insere a nova probabilidade
    Nos = [Nos(1:id-1), {novoNo}, Nos(id:end)]; % Insere o novo nó
end

Arvore = Nos{1}; % Árvore de Huffman
dicionario = mapear(Arvore, ''); % Cria o dicionário de Huffman

TextoCod = codTexto(Dados, dicionario); % Codifica o texto usando Huffman

%% 

% Define matrizes de código de Hamming
H = [ 1 1 0 1 1 0 1 0 1 0 0 0 0;
    1 0 1 1 0 1 1 0 0 1 0 0 0;
    0 1 1 1 0 0 0 1 0 0 1 0 0;
    0 0 0 0 1 1 1 1 0 0 0 1 0;
    1 1 1 0 1 1 0 1 0 0 0 0 1];


c = size(H, 2); % Comprimento do código codificado
d = c - size(H, 1); % Comprimento dos bits de dados

% Define a matriz geradora G
% G =[ I | P' ]
I_k = eye(d);
P = H(:, 1:d); % Extrai a parte da matriz H que representa P
G = [I_k, P'];

% Converte string codificada de Huffman para bits
bits = double(TextoCod) - double('0');
bits = bits(:); % Transpõe para vetor coluna

% Adiciona padding para que o comprimento seja múltiplo de 8 (para o código Hamming (13,8))

% Codifica usando código Hamming (13,8)
dados = reshape(bits, d, [])'; %Faz fatias de 8 bits e a transposição
palavra_codigo = mod(double(dados) * G, 2); % Geração das palavras código
palavra_codigo = palavra_codigo';

% Converte os bits codificados em um vetor coluna
palavra_codigo = palavra_codigo(:);

%%
% Modulação 4-PSK
N = length(palavra_codigo);
cod_mod = zeros(N/2, 1);
for i = 1:2:N % Gera números ímpar    
    par = palavra_codigo(i:i+1); % Faz a fatia de 2 bits. Ex: i=1 -> (1:2), i=3 -> (3:4)
    if isequal(par, [1; 1])
        cod_mod((i+1)/2) = (1 + 1i); % 11 -> (1, 1)
    elseif isequal(par, [1; 0])
        cod_mod((i+1)/2) = (1 - 1i); % 10 -> (1, -1)
    elseif isequal(par, [0; 1])
        cod_mod((i+1)/2) = (-1 + 1i); % 01 -> (-1, 1)
    else
        cod_mod((i+1)/2) = (-1 - 1i); % 00 -> (-1, -1)
    end
end

% Plotar a constelação
figure;
scatter(real(cod_mod), imag(cod_mod), 'filled');
title('Gráfico de Constelação 4-PSK');
xlabel('Parte Real');
ylabel('Parte Imaginária');
xlim([-6 6])
ylim([-6 6])
grid on;

% Demodulação 4-PSK 
recebido = zeros(N, 1);
pontos_constel = [1+1i, 1-1i, -1+1i, -1-1i];
bit_map = [1 1; 1 0; 0 1; 0 0];

% Adicionar ruído e realizar demodulação para diferentes Es/No
Es_No_dB = 0:11; % Valores de Es/No em dB
num_snr = length(Es_No_dB);

ber_sem_correcao = zeros(length(Es_No_dB), 1);
ber_com_correcao = zeros(length(Es_No_dB), 1);
num_erros_matriz = zeros(length(Es_No_dB), 2); % Matriz para registrar número de erros

num_text_erros = zeros(length(Es_No_dB), 1);
for k = 1:num_snr
    % Adicionar ruído usando a função add_awgn_noise
    cod_ruido = add_awgn_noise(cod_mod, Es_No_dB(k));
    
    % Plotar a constelação com ruído
    figure;
    scatter(real(cod_ruido), imag(cod_ruido),1, 'filled');
    hold on;
    scatter(real(pontos_constel), imag(pontos_constel), 5, 'r',  'filled'); 
    title(['Gráfico de Constelação 4-PSK com Ruído (Es/No = ' num2str(Es_No_dB(k)) ' dB)']);
    xlabel('Parte Real');
    ylabel('Parte Imaginária');
    xlim([-6 6]);
    ylim([-6 6]);
    grid on;
    hold off;
    
    % Demodulação 4-PSK
    recebido = zeros(N, 1);
    for i = 1:length(cod_ruido)
        % Encontra o ponto da constelação mais próximo do símbolo recebido
        [~, index] = min(abs(cod_ruido(i) - pontos_constel));
        % Mapeia o índice do ponto encontrado para o par de bits correspondente
        recebido(2*i-1:2*i) = bit_map(index, :);
    end
    
      
    % Calcular BER antes da correção
    uncorrected_bits = reshape(recebido, c, []).';
    uncorrected_bits = uncorrected_bits(:, 1:d).';
    uncorrected_bits = uncorrected_bits(:);
    uncorrected_bits = uncorrected_bits(1:length(bits));
    
    num_errors_no_correction = sum(uncorrected_bits ~= bits);
    ber_sem_correcao(k) = num_errors_no_correction / length(bits);


    % Reshape os bits recebidos para decodificação de Hamming
    entrada_hamming = reshape(recebido, c, []).';
    
    % Calcular o síndrome
    sindrome = mod(entrada_hamming * H.', 2);

    % Corrigir erros com base no síndrome
    for i = 1:size(sindrome, 1)
        if any(sindrome(i, :))
            for j = 1:size(H, 2)
                if isequal(sindrome(i, :), H(:, j)') %Comparacao do sindrome com as colunas de H
                    entrada_hamming(i, j) = mod(entrada_hamming(i, j) + 1, 2); %Flip de correcao
                    break;
                end
            end
        end
    end

    % Extrair bits originais corrigidos
    corrected_bits = entrada_hamming(:, 1:d).';
    corrected_bits = corrected_bits(:);

    % Ajustar o tamanho dos bits corrigidos para corresponder ao tamanho dos bits originais
    corrected_bits = corrected_bits(1:length(bits));

    % Calcular BER após correção
    num_errors_with_correction = sum(corrected_bits ~= bits);
    ber_com_correcao(k) = num_errors_with_correction / length(bits);

    % Registrar o número de erros na matriz
    num_erros_matriz(k, 1) = num_errors_no_correction;
    num_erros_matriz(k, 2) = num_errors_with_correction;

    % Decodificar o texto corrigido usando Huffman
    TextoCod = char(corrected_bits' + '0');
%     TextoDecod = decodTexto(TextoCod, dicionario);

    %Salvar o texto decodificado
%     Decodificado = fopen(['Decodificado_' num2str(Es_No_dB(k)) '.txt'], 'w');
%     disp(['Número de erros no texto decodificado para Es/No = ' num2str(Es_No_dB(k)) ' dB: ' num2str(num_text_errors(k))]);
%     fprintf(Decodificado, '%s', TextoDecod);
%     fclose(Decodificado);

    disp(['Texto decodificado salvo em Decodificado_' num2str(Es_No_dB(k)) '.txt']);
end


% Exibir matriz de erros
disp('Matriz de erros bits (Número de Erros sem Correção, Número de Erros com Correção):');
disp(array2table(num_erros_matriz, 'VariableNames', {'Erros_sem_Correcao', 'Erros_com_Correcao'}, 'RowNames', cellstr(num2str(Es_No_dB(:)))));

figure; % Aumenta a largura da figura
semilogy(Es_No_dB, ber_sem_correcao, '-o', 'Color', [1 0.5 0], 'DisplayName', 'Sem Correção');
xlabel('Es/No (dB)');
ylabel('Taxa de Erro de Bits (BER)');
title('Gráfico de BER sem Correção');
legend('Location', 'best');
grid on;
pbaspect([4 1 1]); 

% Gráfico de BER com Correção
figure; % Aumenta a largura da figura
semilogy(Es_No_dB, ber_com_correcao, '-x', 'DisplayName', 'Com Correção');
xlabel('Es/No (dB)');
ylabel('Taxa de Erro de Bits (BER)');
title('Gráfico de BER com Correção');
legend('Location', 'best');
grid on;
pbaspect([4 1 1]); 
% Gráfico Combinado
figure; % Aumenta a largura da figura
semilogy(Es_No_dB, ber_com_correcao, '-x', 'DisplayName', 'Com Correção');
hold on;
semilogy(Es_No_dB, ber_sem_correcao, '-o', 'DisplayName', 'Sem Correção');
xlabel('Es/No (dB)');
ylabel('Taxa de Erro de Bits (BER)');
title('Gráfico de BER com e sem Correção');
legend('Location', 'best');
grid on;
pbaspect([4 1 1]); 
disp('Processo concluído.');
%% 

% Funções Auxiliares

function checa = JaExiste(Dado, B) 
    n = 1;
    while n <= length(B)
        if Dado == B(n) 
            checa = false;
            return;
        end
        n = n + 1;
    end
    checa = true;
end

function mapa = mapear(arvore, prefixo)
    if iscell(arvore)
        mapa = [mapear(arvore{1}, [prefixo '0']); mapear(arvore{2}, [prefixo '1'])];
    else
        mapa = struct('simbolo', arvore, 'codigo', prefixo);
    end
end

function codificar = codTexto(Dados, dicionario)
    codificar = [];
    for i = 1:length(Dados)
        simbolo = Dados(i);
        id = find([dicionario.simbolo] == simbolo);
        codificar = [codificar dicionario(id).codigo];
    end
end

function texto = decodTexto(codificar, dicionario)
    texto = '';
    codigo = '';
    for i = 1:length(codificar)
        codigo = [codigo codificar(i)]; 
        id = find(strcmp({dicionario.codigo}, codigo));
        if ~isempty(id)
            texto = [texto dicionario(id).simbolo]; 
            codigo = ''; 
        end
    end
end

% Função para adicionar ruído AWGN
function [r, n, N0] = add_awgn_noise(s, SNRdB, L)
    s_temp = s;
    if iscolumn(s), s = s.'; end % Para retornar o resultado no mesmo formato de 's'
    gamma = 10^(SNRdB/10); % SNR para escala linear

    if nargin == 2, L = 1; end % Se o terceiro argumento não for dado, define como 1

    if isvector(s)
        P = L * sum(abs(s).^2) / length(s); % Potência real no vetor
    else % Para sinais multidimensionais como MFSK
        P = L * sum(sum(abs(s).^2)) / length(s); % Se s é uma matriz [MxN]
    end

    N0 = P / gamma; % Encontra a densidade espectral do ruído
    if isreal(s)
        n = sqrt(N0 / 2) * randn(size(s)); % Ruído calculado
    else
        n = sqrt(N0 / 2) * (randn(size(s)) + 1i * randn(size(s))); % Ruído calculado
    end

    r = s + n; % Sinal recebido

    if iscolumn(s_temp), r = r.'; end % Retorna r no formato original como s
end
