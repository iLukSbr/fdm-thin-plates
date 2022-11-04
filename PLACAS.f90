!##########################################################################################
!########     THIN PLATES STRUCTURAL ANALYSIS BY THE FINITE DIFFERENCE METHOD     #########
!##########################################################################################
!###                                                                                    ###
!###    Description: This software has the intent of calculating displacements,         ###
!###    stresses and strains induced in structures.                                     ###
!###                                                                                    ###
!###    Copyright © 2019-2022, Felipe Ferreira de Souza & Lucas Yukio Fukuda Matsumoto  ###
!###    Universidade Federal do Paraná, Centro Politécnico, Setor de Tecnologia,        ###
!###    Departamento de Construção Civil, Curso de Engenharia Civil,                    ###
!###    Av. Cel. Francisco Heráclito dos Santos S/N, Jardim das Américas, Curitiba-PR   ###
!###    Contact by e-mail: lucasyfm@ufpr.br                                             ###
!###                                                                                    ###
!###    This program is free software: you can redistribute it and/or modify            ###
!###    it under the terms of the GNU General Public License as published by            ###
!###    the Free Software Foundation, either version 3 of the License, or               ###
!###    (at your option) any later version.                                             ###
!###                                                                                    ###
!###    This program is distributed in the hope that it will be useful,                 ###
!###    but WITHOUT ANY WARRANTY; without even the implied warranty of                  ###
!###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   ###
!###    GNU General Public License for more details.                                    ###
!###                                                                                    ###
!###    You should have received a copy of the GNU General Public License               ###
!###    along with this program.  If not, see <https://www.gnu.org/licenses/>           ###
!###                                                                                    ###
!##########################################################################################

PROGRAM PLACAS_DELGADAS

IMPLICIT NONE

!VARIÁVEIS
INTEGER :: a, b, i, j, k, l, m, n, qtd_bordas, x_div, y_div, num_div, ierror
INTEGER :: n_zeros, linhas_nulas, colunas_nulas, pontos, pontos_extra, info, valor, pt_i, pt_jme1
INTEGER :: pt_j, pt_jma1, pt_kme2, pt_kme1, pt_k, pt_kma1, pt_kma2, pt_lme1, pt_l, pt_lma1, pt_m
REAL :: comprimento, largura, espessura, menordim, elast_long, delta_x, delta_y, inercia_x, inercia_y, erro_relativo
REAL :: poisson, carga, mod_rigidez, wmax_mdf, xmax_mdf, ymax_mdf, wmax_calc, xmax_calc, ymax_calc, somatorio, pi, aux
REAL :: mxmax_calc, mymax_calc, beta_m, A_m, B_m, C_m, D_m
CHARACTER(1) :: opcao
CHARACTER(2) :: tipo_grafico, tipo_arquivo, direcao
CHARACTER(4), DIMENSION(:), ALLOCATABLE :: tipo_borda, posicao_borda
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
REAL, DIMENSION(:), ALLOCATABLE :: vetor_desloc, vetor_abscissa, vetor_ordenada, vet_cargas
REAL, DIMENSION(:,:), ALLOCATABLE :: mat_coeficientes, mat_pontos, mat_deslocamento, mat_expandida, mat_fx, mat_fy

909 CONTINUE
WRITE(*,*) "####################################################################################"
WRITE(*,*) "####  ANÁLISE ESTRUTURAL DE PLACAS DELGADAS PELO MÉTODO DAS DIFERENÇAS FINITAS  ####"
WRITE(*,*) "####################################################################################"
WRITE(*,*) " "
WRITE(*,*) "                       ______________________________________"
WRITE(*,*) "                      /                                     /|"
WRITE(*,*) "                     /                                     / /"
WRITE(*,*) "                    /                                     / /"
WRITE(*,*) "                   /           PLACA DELGADA             / /"
WRITE(*,*) "       LARGURA  Y /              SUJEITA A              / /"
WRITE(*,*) "                 /         PEQUENAS DEFLEXÕES          / /"
WRITE(*,*) "                /                                     / /"
WRITE(*,*) "               /                                     / /"
WRITE(*,*) "              /                                     / /"
WRITE(*,*) "             /_____________________________________/ /"
WRITE(*,*) "ESPESSURA  Z |_____________________________________|/"
WRITE(*,*) "                                X"
WRITE(*,*) "                           COMPRIMENTO"
WRITE(*,*) " "
WRITE(*,*) "> Siga as instruções e solicitações;"
WRITE(*,*) "> Utilize o ponto (.) como separador decimal;"
WRITE(*,*) "> Não utilize separador de milhar;"
WRITE(*,*) "> Não digite as unidades de medida;"
WRITE(*,*) "> Não insira espaços;"
WRITE(*,*) "> A tecla Enter conclui a sua digitação."
WRITE(*,*) " "

WRITE(*,*) "teste:PROPRIEDADES GEOMÉTRICAS"
333 CONTINUE
WRITE(*,*) "Digite o comprimento da placa em metros (direção X)."
READ(*,*) comprimento
IF(comprimento <= 0) THEN
    WRITE(*,*) "Valor inválido, digite novamente."
    GO TO 333
END IF
334 CONTINUE
WRITE(*,*) "Digite a largura da placa em metros (direção Y)."
READ(*,*) largura
IF(largura <= 0) THEN
    WRITE(*,*) "Valor inválido, digite novamente."
    GO TO 334
END IF
335 CONTINUE
WRITE(*,*) "Digite a espessura da placa em centímetros (direção Z)."
READ(*,*) espessura
IF(espessura <= 0) THEN
    WRITE(*,*) "Valor inválido, digite novamente."
    GO TO 335
END IF
336 CONTINUE
WRITE(*,*) "Digite o número de divisões nas duas direções."
READ(*,*) x_div
IF(x_div <= 3) THEN
    WRITE(*,*) "Escassez de divisões, tente uma precisão mais alta."
    GO TO 336
END IF
y_div = x_div
espessura = 0.01*espessura
qtd_bordas = 4
delta_x = comprimento/(x_div-1)
delta_y = largura/(y_div-1)
inercia_x = largura*(espessura**3)/12
inercia_y = comprimento*(espessura**3)/12

WRITE(*,*) "teste:VERIFICAÇÃO SE A PLACA É DELGADA"
IF(comprimento >= largura) THEN
    menordim = largura
ELSE
    menordim = comprimento
END IF
WRITE(*,*) "A relação 20 < Menor dimensão/Espessura =",menordim/espessura," < 80 deve ser verdadeira para uma placa delgada,"
IF(menordim/espessura > 20 .AND. menordim/espessura < 80) THEN
    WRITE(*,*) "A geometria da placa atende ao critério de placas delgadas."
ELSE IF(menordim/espessura >= 80) THEN
    WRITE(*,*) "A geometria da placa não atende ao critério de placas delgadas, é uma membrana."
    GO TO 333
ELSE IF(menordim/espessura <= 20) THEN
    WRITE(*,*) "A geometria da placa não atende ao critério de placas delgadas, é uma placa espessa."
    GO TO 333
END IF

WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO"
ALLOCATE(tipo_borda(qtd_bordas))
ALLOCATE(posicao_borda(qtd_bordas))
WRITE(*,*) " "
WRITE(*,*) "BSA = Borda Simplesmente Apoiada"
WRITE(*,*) "BE = Borda Engastada"
WRITE(*,*) "BL = Borda Livre"
WRITE(*,*) "BD = Borda Deslizante"
WRITE(*,*) " "
WRITE(*,*) "> Seta para cima repete a última digitação."
WRITE(*,*) "> Para voltar ao início, digite 'V'."
WRITE(*,*) "> Para sair, digite 'S'."
WRITE(*,*) " "
710 CONTINUE
i = 1
WRITE(*,*) "Qual a condição da borda superior?"
READ(*,*) tipo_borda(i)
IF(tipo_borda(i) == 'V') THEN
    DEALLOCATE(tipo_borda)
    DEALLOCATE(posicao_borda)
    GO TO 909
ELSE IF(tipo_borda(i) == 'S') THEN
    STOP
ELSE IF(tipo_borda(i) /= 'BSA' .AND. tipo_borda(i) /= 'BE' .AND. tipo_borda(i) /= 'BL' .AND. tipo_borda(i) /= 'BD') THEN
    WRITE(*,*) "Tipo inexistente."
    GO TO 710
END IF
posicao_borda(i) = 'SUP'
i = i+1
720 CONTINUE
WRITE(*,*) "Qual a condição da borda direita?"
READ(*,*) tipo_borda(i)
IF(tipo_borda(i) == 'V') THEN
    DEALLOCATE(tipo_borda)
    DEALLOCATE(posicao_borda)
    GO TO 909
ELSE IF(tipo_borda(i) == 'S') THEN
    STOP
ELSE IF(tipo_borda(i) /= 'BSA' .AND. tipo_borda(i) /= 'BE' .AND. tipo_borda(i) /= 'BL' .AND. tipo_borda(i) /= 'BD') THEN
    WRITE(*,*) "Tipo inexistente."
    GO TO 720
END IF
posicao_borda(i) = 'DIR'
i = i+1
730 CONTINUE
WRITE(*,*) "Qual a condição da borda inferior?"
READ(*,*) tipo_borda(i)
IF(tipo_borda(i) == 'V') THEN
    DEALLOCATE(tipo_borda)
    DEALLOCATE(posicao_borda)
    GO TO 909
ELSE IF(tipo_borda(i) == 'S') THEN
    STOP
ELSE IF(tipo_borda(i) /= 'BSA' .AND. tipo_borda(i) /= 'BE' .AND. tipo_borda(i) /= 'BL' .AND. tipo_borda(i) /= 'BD') THEN
    WRITE(*,*) "Tipo inexistente."
    GO TO 730
END IF
posicao_borda(i) = 'INF'
i = i+1
740 CONTINUE
WRITE(*,*) "Qual a condição da borda esquerda?"
READ(*,*) tipo_borda(i)
IF(tipo_borda(i) == 'V') THEN
    DEALLOCATE(tipo_borda)
    DEALLOCATE(posicao_borda)
    GO TO 909
ELSE IF(tipo_borda(i) == 'S') THEN
    STOP
ELSE IF(tipo_borda(i) /= 'BSA' .AND. tipo_borda(i) /= 'BE' .AND. tipo_borda(i) /= 'BL' .AND. tipo_borda(i) /= 'BD') THEN
    WRITE(*,*) "Tipo inexistente."
    GO TO 740
END IF
posicao_borda(i) = 'ESQ'
j = 0
k = 0
m = 0
DO i = 1,qtd_bordas
    IF(tipo_borda(i) == 'BSA') THEN
        j = j+1
    ELSE IF(tipo_borda(i) == 'BD') THEN
        k = k+1
    ELSE IF(tipo_borda(i) == 'BL') THEN
        m = m+1
    END IF
END DO
IF(k == 4 .OR. m == 4 .OR. (j == 1 .AND. m == 3)) THEN
    WRITE(*,*) "Placa sem sustentação (4xBD ou 4xBL ou 1xBSA + 3xBL), corrigir as bordas."
    GO TO 710
END IF
DO i = 1,qtd_bordas
    IF(tipo_borda(i) == 'BSA' .AND. (posicao_borda(i) == 'SUP' .OR. posicao_borda(i) == 'INF')) THEN
        tipo_borda(i) = 'BSAx'
    ELSE IF(tipo_borda(i) == 'BSA' .AND. (posicao_borda(i) == 'ESQ' .OR. posicao_borda(i) == 'DIR')) THEN
        tipo_borda(i) = 'BSAy'
    ELSE IF(tipo_borda(i) == 'BE' .AND. (posicao_borda(i) == 'SUP' .OR. posicao_borda(i) == 'INF')) THEN
        tipo_borda(i) = 'BEx'
    ELSE IF(tipo_borda(i) == 'BE' .AND. (posicao_borda(i) == 'ESQ' .OR. posicao_borda(i) == 'DIR')) THEN
        tipo_borda(i) = 'BEy'
    ELSE IF(tipo_borda(i) == 'BL' .AND. (posicao_borda(i) == 'SUP' .OR. posicao_borda(i) == 'INF')) THEN
        tipo_borda(i) = 'BLx'
    ELSE IF(tipo_borda(i) == 'BL' .AND. (posicao_borda(i) == 'ESQ' .OR. posicao_borda(i) == 'DIR')) THEN
        tipo_borda(i) = 'BLy'
    ELSE IF(tipo_borda(i) == 'BD' .AND. (posicao_borda(i) == 'SUP' .OR. posicao_borda(i) == 'INF')) THEN
        tipo_borda(i) = 'BDx'
    ELSE IF(tipo_borda(i) == 'BD' .AND. (posicao_borda(i) == 'ESQ' .OR. posicao_borda(i) == 'DIR')) THEN
        tipo_borda(i) = 'BDy'
    END IF
END DO

WRITE(*,*) "teste:PROPRIEDADES DO MATERIAL"
337 CONTINUE
WRITE(*,*) "Digite o módulo de elasticidade longitudinal do material em GPa"
READ(*,*) elast_long
IF(elast_long <= 0) THEN
    WRITE(*,*) "Valor inválido, digite novamente."
    GO TO 337
END IF
elast_long = 1000000000*elast_long
WRITE(*,*) "Digite o coeficiente de Poisson do material"
READ(*,*) poisson
338 CONTINUE
IF(poisson <= 0) THEN
    WRITE(*,*) "Valor inválido, digite novamente."
    GO TO 338
END IF
mod_rigidez = elast_long*(espessura**3)/(12*(1 - poisson**2))

WRITE(*,*) "teste:CARREGAMENTO"
339 CONTINUE
WRITE(*,*) "Digite o valor da carga distribuída em kN/m²"
READ(*,*) carga
IF(carga <= 0) THEN
    WRITE(*,*) "Valor inválido, digite novamente."
    GO TO 339
END IF
carga = 1000*carga
num_div = x_div*y_div

WRITE(*,*) "teste:TAMANHO E ALOCAÇÃO DA MATRIZ DE PONTOS E REMOÇÃO DE RESÍDUOS"
pontos = x_div*y_div
pontos_extra = (x_div+4)*(y_div+4)
WRITE(*,*) "teste:Tamanho da matriz de pontos com os extras:",pontos,"x",pontos_extra
ALLOCATE(mat_pontos(pontos,pontos_extra))
DO i = 1, SIZE(mat_pontos,1)
    DO j = 1, SIZE(mat_pontos,2)
        mat_pontos(i,j) = 0
    END DO
END DO

WRITE(*,*) "teste:APLICAÇÃO DA EQUAÇÃO GOVERNANTE EM TODOS OS PONTOS"
j = 2*x_div+11
b = 1
a = 1
DO i = 1, pontos
! WRITE(*,*) "teste:Ponto na linha",i,"e coluna",j
    CALL GERMAIN_LAGRANGE
    j = j+1
    a = a+1
    IF(a > x_div) THEN
        j = 2*x_div+11 + b*(x_div+4) !Colunas andadas (pontos) no total
        b = b+1 !Linhas andadas (pontos)
        a = 1 !Colunas andadas (pontos) na mesma linha
    END IF
END DO
OPEN(UNIT=10, FILE='GermainLagrange.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, SIZE(mat_pontos,1)
        WRITE(10,*) mat_pontos(i,:)
    END DO
CLOSE (10)

DO k = 1, qtd_bordas
    WRITE(*,*) "teste:APLICA AS CONDIÇÕES DE BORDA ",tipo_borda(k)," ",posicao_borda(k)
    IF(tipo_borda(k) == 'BSAx' .OR. tipo_borda(k) == 'BSAy') THEN
        CALL BSA
    ELSE IF(tipo_borda(k) == 'BEx' .OR. tipo_borda(k) == 'BEy') THEN
        CALL BE
    ELSE IF(tipo_borda(k) == 'BDx' .OR. tipo_borda(k) == 'BDy') THEN
        CALL BD
        CALL BE
    ELSE IF(tipo_borda(k) == 'BLx' .OR. tipo_borda(k) == 'BLy') THEN
        CALL BL
        CALL BSA
    END IF
END DO

DO k = 1, qtd_bordas
    IF(k < qtd_bordas) THEN
        a = k+1
    ELSE
        a = 1
    END IF
    WRITE(*,*) "teste:APLICA AS CONDIÇÕES DE CANTO ",tipo_borda(k)," ",tipo_borda(a)
    IF((tipo_borda(k) == 'BSAx' .AND. tipo_borda(a) =='BSAy') .OR. &
       (tipo_borda(k) == 'BSAy' .AND. tipo_borda(a) == 'BSAx') .OR. &
       (tipo_borda(k) == 'BEx' .AND. tipo_borda(a) == 'BEy') .OR. &
       (tipo_borda(k) == 'BEy' .AND. tipo_borda(a) == 'BEx') .OR. &
       (tipo_borda(k) == 'BSAx' .AND. tipo_borda(a) == 'BEy') .OR. &
       (tipo_borda(k) == 'BEx' .AND. tipo_borda(a) == 'BSAy') .OR. &
       (tipo_borda(k) == 'BSAy' .AND. tipo_borda(a) == 'BEx') .OR. &
       (tipo_borda(k) == 'BEy' .AND. tipo_borda(a) == 'BSAx')) THEN
        CALL CANTO_BSABSA_BEBE_BSABE
    ELSE IF((tipo_borda(k) == 'BSAx' .AND. tipo_borda(a) == 'BDy') .OR. &
            (tipo_borda(k) == 'BDx' .AND. tipo_borda(a) == 'BSAy') .OR. &
            (tipo_borda(k) == 'BSAy' .AND. tipo_borda(a) == 'BDx') .OR. &
            (tipo_borda(k) == 'BDy' .AND. tipo_borda(a) == 'BSAx')) THEN
        CALL CANTO_BSABD
    ELSE IF((tipo_borda(k) == 'BDx' .AND. tipo_borda(a) == 'BLy') .OR. &
            (tipo_borda(k) == 'BLx' .AND. tipo_borda(a) == 'BDy') .OR. &
            (tipo_borda(k) == 'BDy' .AND. tipo_borda(a) == 'BLx') .OR. &
            (tipo_borda(k) == 'BLy' .AND. tipo_borda(a) == 'BDx')) THEN
        CALL CANTO_BDBL
    ELSE IF((tipo_borda(k) == 'BSAx' .AND. tipo_borda(a) == 'BLy') .OR. &
            (tipo_borda(k) == 'BLx' .AND. tipo_borda(a) == 'BSAy') .OR. &
            (tipo_borda(k) == 'BSAy' .AND. tipo_borda(a) == 'BLx') .OR. &
            (tipo_borda(k) == 'BLy' .AND. tipo_borda(a) == 'BSAx')) THEN
        CALL CANTO_BSABL
    ELSE IF((tipo_borda(k) == 'BSAx' .AND. tipo_borda(a) == 'BLy') .OR. &
            (tipo_borda(k) == 'BLx' .AND. tipo_borda(a) == 'BSAy') .OR. &
            (tipo_borda(k) == 'BSAy' .AND. tipo_borda(a) == 'BLx') .OR. &
            (tipo_borda(k) == 'BLy' .AND. tipo_borda(a) == 'BSAx')) THEN
        CALL CANTO_BSABL
    ELSE IF((tipo_borda(k) == 'BDx' .AND. tipo_borda(a) == 'BDy') .OR. &
            (tipo_borda(k) == 'BDy' .AND. tipo_borda(a) == 'BDx')) THEN
        CALL CANTO_BDBD
    ELSE IF((tipo_borda(k) == 'BLx' .AND. tipo_borda(a) == 'BLy') .OR. &
            (tipo_borda(k) == 'BLy' .AND. tipo_borda(a) == 'BLx')) THEN
        CALL CANTO_BLBL
    END IF
END DO
OPEN(UNIT=10, FILE='CoeficienteComNulos.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, SIZE(mat_pontos,1)
        WRITE(10,*) mat_pontos(i,:)
    END DO
CLOSE (10)

WRITE(*,*) "teste:CONTA LINHAS ZERADAS"
linhas_nulas = 0
DO i = 1, pontos
    n_zeros = 0
    DO j = 1, pontos_extra
        IF(mat_pontos(i,j) == 0) THEN
            n_zeros = n_zeros+1
        END IF
    END DO
    IF(n_zeros == pontos_extra) THEN
        linhas_nulas = linhas_nulas+1
    END IF
END DO

WRITE(*,*) "teste:CONTA COLUNAS ZERADAS"
colunas_nulas = 0
DO j = 1, pontos_extra
    n_zeros = 0
    DO i = 1, pontos
        IF(mat_pontos(i,j) == 0) THEN
            n_zeros = n_zeros+1
        END IF
    END DO
    IF(n_zeros == pontos) THEN
        colunas_nulas = colunas_nulas+1
    END IF
END DO

WRITE(*,*) "teste:Tamanho da matriz de coeficientes:",pontos-linhas_nulas,"x",pontos_extra-colunas_nulas
WRITE(*,*) "teste:ALOCA E ELIMINA RESÍDUOS DA MATRIZ DE COEFICIENTES SEM EQUAÇÕES NULAS"
ALLOCATE(mat_coeficientes(pontos-linhas_nulas, pontos_extra-colunas_nulas))
DO i = 1, SIZE(mat_coeficientes,1)
    DO j = 1, SIZE(mat_coeficientes,2)
        mat_coeficientes(i,j) = 0
    END DO
END DO

WRITE(*,*) "teste:ARMAZENA A MATRIZ DE COEFICIENTES APARADA"
k = 1
DO i = 1, pontos
    valor = 0
    l = 1
    DO j = 1, pontos_extra
        n_zeros = 0
        DO m = 1, pontos
            IF(mat_pontos(m,j) == 0) THEN
                n_zeros = n_zeros+1
            END IF
            IF(n_zeros < pontos) THEN
                a = 1 !A coluna não é nula
            ELSE
                a = 0 !A coluna é nula
            END IF
        END DO
        n_zeros = 0
        DO n = 1, pontos_extra
            IF(mat_pontos(i,n) == 0) THEN
                n_zeros = n_zeros+1
            END IF
            IF(n_zeros < pontos_extra) THEN
                b = 1 !A linha não é nula
            ELSE
                b = 0 !A linha é nula
            END IF
        END DO
        IF(a == 1 .AND. b == 1) THEN
            valor = valor+1 !Teve valor copiado nessa linha
            mat_coeficientes(k,l) = mat_pontos(i,j)
            l = l+1
        END IF
    END DO
    IF(valor > 0) THEN
        k = k+1
        valor = 0
    END IF
END DO
OPEN(UNIT=10, FILE='Coeficientes.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, SIZE(mat_coeficientes,1)
        WRITE(10,*) mat_coeficientes(i,:)
    END DO
CLOSE (10)

WRITE(*,*) "teste:RESOLUÇÃO DO SISTEMA LINEAR"
m = SIZE(mat_coeficientes,1)
ALLOCATE(ipiv(m))
ALLOCATE(vet_cargas(m))
DO j = 1, SIZE(vet_cargas,1) !vetor Pk/D
    vet_cargas(j) = carga/mod_rigidez
END DO
CALL DGESV(m, 1, mat_coeficientes, m, ipiv, vet_cargas, m, info)
vetor_desloc = vet_cargas
OPEN(UNIT=10, FILE='Vetor Deslocamento.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, SIZE(vetor_desloc,1)
        WRITE(10,*) vetor_desloc(i)
    END DO
CLOSE (10)

WRITE(*,*) "teste:ELABORA LISTAS DE COORDENADAS X E Y DOS PONTOS"
ALLOCATE(vetor_abscissa(num_div))
ALLOCATE(vetor_ordenada(num_div))
k = 1
DO i = 1, y_div
    DO j = 1, x_div
        vetor_abscissa(k) = (j-1)*delta_x
        vetor_ordenada(k) = (i-1)*delta_y
        k = k+1
    END DO
END DO

WRITE(*,*) "teste:ALOCA A MATRIZ DE DESLOCAMENTOS, ELIMINA RESÍDUOS E PREENCHE PARA CADA PONTO DA PLACA"
ALLOCATE (mat_deslocamento(y_div,x_div))
DO i = 1, SIZE(mat_deslocamento,1)
    DO j = 1, SIZE(mat_deslocamento,2)
        mat_deslocamento(i,j) = 0
    END DO
END DO
k = 1
DO i = 1, y_div
    IF((tipo_borda(1) == 'BSAx' .OR. tipo_borda(1) == 'BEx') .AND. &
       (tipo_borda(4) == 'BSAy' .OR. tipo_borda(4) == 'BEy')) THEN
        IF((tipo_borda(1) == 'BSAx' .OR. tipo_borda(1) == 'BEx') .AND. &
           (tipo_borda(2) == 'BSAy' .OR. tipo_borda(2) == 'BEy')) THEN
            DO j = 2, x_div-1 !Borda superior, esquerda e direita w=0 (000)
                mat_deslocamento(i+1,j) = vetor_desloc(k)
                k = k+1
            END DO
        ELSE
            DO j = 2, x_div !Borda superior e esquerda w=0 e Borda direita w>0 (001)
                mat_deslocamento(i+1,j) = vetor_desloc(k)
                k = k+1
            END DO
        END IF
    ELSE IF((tipo_borda(1) == 'BSAx' .OR. tipo_borda(1) == 'BEx') .AND. &
           (tipo_borda(2) == 'BSAy' .OR. tipo_borda(2) == 'BEy')) THEN
        DO j = 1, x_div-1 !Borda esquerda w>0 e Borda superior e direita w=0 (100)
            mat_deslocamento(i+1,j) = vetor_desloc(k)
            k = k+1
        END DO
    ELSE IF(tipo_borda(4) == 'BSAy' .OR. tipo_borda(4) == 'BEy') THEN
        IF(tipo_borda(2) == 'BSAy' .OR. tipo_borda(2) == 'BEy') THEN
            DO j = 2, x_div-1 !Borda esquerda e direita w=0 e Borda superior w>0 (010)
                mat_deslocamento(i,j) = vetor_desloc(k)
                k = k+1
            END DO
        ELSE
            DO j = 2, x_div !Borda esquerda w=0 e Borda superior e direita w>0 (011)
                mat_deslocamento(i,j) = vetor_desloc(k)
                k = k+1
            END DO
        END IF
    ELSE IF(tipo_borda(1) == 'BSAx' .OR. tipo_borda(1) == 'BEx') THEN
        DO j = 1, x_div !Borda esquerda e direita w>0 e Borda superior w=0 (101)
            mat_deslocamento(i+1,j) = vetor_desloc(k)
            k = k+1
        END DO
    ELSE IF(tipo_borda(2) == 'BSAy' .OR. tipo_borda(2) == 'BEy') THEN
        DO j = 1, x_div-1 !Borda esquerda e superior w>0 e Borda direita w=0 (110)
            mat_deslocamento(i,j) = vetor_desloc(k)
            k = k+1
        END DO
    ELSE !Borda superior, esquerda e direita w>0 (111)
        DO j = 1, x_div
            mat_deslocamento(i,j) = vetor_desloc(k)
            k = k+1
        END DO
    END IF
    IF((k > SIZE(vetor_desloc,1)) .OR. ((tipo_borda(3) == 'BSAx' .OR. tipo_borda(3) == 'BEx') .AND. &
       i == y_div-1) .OR. ((tipo_borda(1) == 'BSAx' .OR. tipo_borda(1) == 'BEx') .AND. & 
       (tipo_borda(3) == 'BSAx' .OR. tipo_borda(3) == 'BEx') .AND. i == y_div-2)) THEN
        EXIT
    END IF
END DO
k = 1
OPEN(UNIT=10, FILE='Deslocamento.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, y_div
        DO j = 1, x_div
            WRITE(10,*) vetor_abscissa(k), vetor_ordenada(k), 1000*mat_deslocamento(i,j)
            k = k+1
        END DO
        WRITE(10,*) " "
    END DO
CLOSE (10)
DO i = 1, SIZE(mat_deslocamento,1)
    DO j = 1, SIZE(mat_deslocamento,2)
        IF(mat_deslocamento(i,j) > wmax_mdf) THEN
            wmax_mdf = mat_deslocamento(i,j)
            xmax_mdf = (j-1)*(comprimento/(x_div-1))
            ymax_mdf = (i-1)*(largura/(y_div-1))
        END IF
    END DO
END DO
WRITE(*,*) "Deslocamento máximo por MDF de ",wmax_mdf*1000,"mm em X=",xmax_mdf," Y=",ymax_mdf
IF(wmax_mdf <= 0.2*espessura) THEN
    WRITE(*,*) "Valor válido para a teoria de pequenas deflexões, o deslocamento é pequeno:"
    WRITE(*,*) "w=",wmax_mdf*1000,"mm <= 0,2h=",0.2*espessura*1000,"mm"
ELSE
    WRITE(*,*) "Valor INVÁLIDO para a teoria de pequenas deflexões, o deslocamento é grande:"
    WRITE(*,*) "w=",wmax_mdf*1000,"mm > 0,2h=",0.2*espessura*1000,"mm"
    WRITE(*,*) "Digite 'V' para voltar ao início ou 'S' para sair."
    READ(*,*) opcao
    IF(opcao == 'V') THEN
    ELSE IF(tipo_arquivo == 'V') THEN
        DEALLOCATE(tipo_borda)
        DEALLOCATE(posicao_borda)
        DEALLOCATE(mat_pontos)
        DEALLOCATE(mat_coeficientes)
        DEALLOCATE(ipiv)
        DEALLOCATE(vet_cargas)
        DEALLOCATE(vetor_ordenada)
        DEALLOCATE(vetor_abscissa)
        DEALLOCATE(mat_deslocamento)
        GO TO 909
    ELSE
        STOP
    END IF
END IF

WRITE(*,*) "teste:ALOCA E PREENCHE A MATRIZ DE DESLOCAMENTOS EXPANDIDA"
ALLOCATE(mat_expandida(SIZE(mat_deslocamento,1)+4,SIZE(mat_deslocamento,2)+4))
DO i = 1, SIZE(mat_expandida,1)
    DO j = 1, SIZE(mat_expandida,2)
        mat_expandida(i,j) = 0
    END DO
END DO
m = 1
DO i = 3, SIZE(mat_expandida,1)-2
    n = 1
    DO j = 3, SIZE(mat_expandida,2)-2
        mat_expandida(i,j) = mat_deslocamento(m,n)
        n = n+1
    END DO
    m = m+1
END DO
DO k = 1, qtd_bordas
    CALL EXPANSAO
END DO
ALLOCATE(mat_fx(x_div,y_div))
CALL F_x
ALLOCATE(mat_fy(x_div,y_div))
CALL F_y
k = 1
OPEN(UNIT=20, FILE='Momento Fletor X.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, y_div
        DO j = 1, x_div
            WRITE(20,*) vetor_abscissa(k), vetor_ordenada(k), mat_fx(i,j)/1000
            k = k+1
        END DO
        WRITE(10,*) " "
    END DO
CLOSE (20)
k = 1
OPEN(UNIT=30, FILE='Momento Fletor Y.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    DO i = 1, y_div
        DO j = 1, x_div
            WRITE(30,*) vetor_abscissa(k), vetor_ordenada(k), mat_fy(i,j)/1000
            k = k+1
        END DO
        WRITE(30,*) " "
    END DO
CLOSE (30)

k = 0
DO i = 1, qtd_bordas
    IF(tipo_borda(i) == 'BSAx' .OR. tipo_borda(i) == 'BSAy') THEN
        k = k+1
    END IF
END DO
IF(k == 4) THEN
    WRITE(*,*) "teste:SOLUÇÃO DE NAVIER"
    somatorio = 0
    pi = 4*ATAN(1.d0)
    xmax_calc = 0.5*comprimento
    ymax_calc = 0.5*largura
    DO m = 1, 5000
        DO n = 1, 5000
            aux = 0.5*(m+n)
            IF(aux-FLOOR(aux) == 0) THEN
                somatorio = somatorio + ((-1)**(aux-1))/(m*n*((m/comprimento)**2+(n/largura)**2)**2)
                ! WRITE(*,*)"teste:",m,n,somatorio
            END IF
            mxmax_calc = mxmax_calc + (((m/comprimento)**2 + poisson*(n/largura)**2)/(m*n*((m/comprimento)**2 &
                                    + (n/largura)**2)**2))*SIN(m*pi*xmax_calc/comprimento)*SIN(n*pi*ymax_calc/largura)
            mymax_calc = mymax_calc + ((poisson*(m/comprimento)**2 + (n/largura)**2)/(m*n*((m/comprimento)**2 &
                                    + (n/largura)**2)**2))*SIN(m*pi*xmax_calc/comprimento)*SIN(n*pi*ymax_calc/largura)
        END DO
    END DO
    wmax_calc = (16*carga/((pi**6)*mod_rigidez))*somatorio
    WRITE(*,*) "Deslocamento máximo por Navier de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
    erro_relativo = ABS(wmax_mdf-wmax_calc)/wmax_calc
    WRITE(*,*) "Erro relativo =",erro_relativo*100,"%"
    mxmax_calc = mxmax_calc*16*carga/(pi**4)
    mymax_calc = mymax_calc*16*carga/(pi**4)
    WRITE(*,*) "Momento Fletor X máximo por Navier de ",mxmax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
    WRITE(*,*) "Momento Fletor Y máximo por Navier de ",mymax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
    OPEN(UNIT=40, FILE='Deslocamento Maximo.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
        WRITE(40,*) "Deslocamento máximo por MDF de ",wmax_mdf*1000,"mm em X=",xmax_mdf," Y=",ymax_mdf
        WRITE(40,*) "Deslocamento máximo por Navier de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
        WRITE(40,*) "Erro relativo =",erro_relativo*100,"%"
        WRITE(40,*) "Momento Fletor X máximo por Navier de ",mxmax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
        WRITE(40,*) "Momento Fletor Y máximo por Navier de ",mymax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
    CLOSE (40)
END IF

IF(tipo_borda(1) == 'BEx' .AND. largura /= comprimento) THEN
    k = 0
    DO i = 2, qtd_bordas
        IF(tipo_borda(i) == 'BSAx' .OR. tipo_borda(i) == 'BSAy') THEN
            k = k+1
        END IF
    END DO
    IF(k == 3) THEN
        WRITE(*,*) "teste:SOLUÇÃO DE LÉVY"
        pi = 4*ATAN(1.d0)
        xmax_calc = 0.5*comprimento
        ymax_calc = ymax_mdf
        DO m = 1, 79, 2
            beta_m = m*pi*largura/comprimento
            A_m = (2*carga*(comprimento**4)/((m**5)*(pi**5)*mod_rigidez))*(2*(COSH(beta_m)**2) - 2*COSH(beta_m) &
                - beta_m*SINH(beta_m))/(COSH(beta_m)*SINH(beta_m)-beta_m)
            B_m = -4*carga*(comprimento**4)/((m**5)*(pi**5)*mod_rigidez)
            C_m = (2*carga*(comprimento**4)/((m**5)*(pi**5)*mod_rigidez))*(beta_m/largura)*(2*SINH(beta_m)*COSH(beta_m) &
                - SINH(beta_m) - beta_m*COSH(beta_m))/(SINH(beta_m)*COSH(beta_m)-beta_m)
            D_m = -A_m*m*pi/comprimento
            wmax_calc = wmax_calc + (A_m*SINH(m*pi*ymax_calc/comprimento) + B_m*COSH(m*pi*ymax_calc/comprimento) + C_m*ymax_calc&
                      *SINH(m*pi*ymax_calc/comprimento) + D_m*ymax_calc*SINH(m*pi*ymax_calc/comprimento) + 4*carga*(comprimento**4)&
                      /((m**5)*(pi**5)*mod_rigidez))*SIN(m*pi*xmax_calc/comprimento)
            WRITE(*,*) "teste:",m,A_m,B_m,C_m,D_m,wmax_calc
        END DO
        WRITE(*,*) "Deslocamento máximo por Lévy de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
        erro_relativo = ABS(wmax_mdf-wmax_calc)/wmax_calc
        WRITE(*,*) "Erro relativo =",erro_relativo*100,"%"
        OPEN(UNIT=40, FILE='Deslocamento Maximo.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
            WRITE(40,*) "Deslocamento máximo por MDF de ",wmax_mdf*1000,"mm em X=",xmax_mdf," Y=",ymax_mdf
            WRITE(40,*) "Deslocamento máximo por Lévy de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
            WRITE(40,*) "Erro relativo =",erro_relativo*100,"%"
        CLOSE (40)
    END IF
END IF

IF(tipo_borda(1) == 'BEx' .AND. largura == comprimento) THEN
    k = 0
    DO i = 2, qtd_bordas
        IF(tipo_borda(i) == 'BSAx' .OR. tipo_borda(i) == 'BSAy') THEN
            k = k+1
        END IF
    END DO
    IF(k == 3) THEN
        WRITE(*,*) "teste:SOLUÇÃO DE LÉVY QUADRADA"
        pi = 4*ATAN(1.d0)
        xmax_calc = 0.5*comprimento
        ymax_calc = ymax_mdf
        wmax_calc = 0.0028*carga*(comprimento**4)/mod_rigidez
        mymax_calc = 0.084*carga*comprimento**2
        WRITE(40,*) "Momento Fletor Y máximo por Superposição de ",mymax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
        WRITE(*,*) "Deslocamento máximo por Lévy de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
        erro_relativo = ABS(wmax_mdf-wmax_calc)/wmax_calc
        WRITE(*,*) "Erro relativo =",erro_relativo*100,"%"
        OPEN(UNIT=40, FILE='Deslocamento Maximo.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
            WRITE(40,*) "Deslocamento máximo por MDF de ",wmax_mdf*1000,"mm em X=",xmax_mdf," Y=",ymax_mdf
            WRITE(40,*) "Deslocamento máximo por Lévy de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
            WRITE(40,*) "Erro relativo =",erro_relativo*100,"%"
        CLOSE (40)
    END IF
END IF

IF(largura == 2*comprimento) THEN
    k = 0
    DO i = 1, qtd_bordas
        IF(tipo_borda(i) == 'BEx' .OR. tipo_borda(i) == 'BEy') THEN
            k = k+1
        END IF
    END DO
    IF(k == 4) THEN
        WRITE(*,*) "teste:SOLUÇÃO PELO MÉTODO DA SUPERPOSIÇÃO"
        xmax_calc = 0.5*comprimento
        ymax_calc = 0.5*largura
        wmax_calc = 0.00254*carga*(comprimento**4)/mod_rigidez
        WRITE(*,*) "Deslocamento máximo por Superoposição de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
        erro_relativo = ABS(wmax_mdf-wmax_calc)/wmax_calc
        WRITE(*,*) "Erro relativo =",erro_relativo*100,"%"
        mxmax_calc = 0.0829*carga*comprimento**2
        mymax_calc = mxmax_calc
        WRITE(40,*) "Momento Fletor X máximo por Superposição de ",mxmax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
        WRITE(40,*) "Momento Fletor Y máximo por Superposição de ",mymax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
        OPEN(UNIT=40, FILE='Deslocamento Maximo.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
            WRITE(40,*) "Deslocamento máximo por MDF de ",wmax_mdf*1000,"mm em X=",xmax_mdf," Y=",ymax_mdf
            WRITE(40,*) "Deslocamento máximo por Superposição de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
            WRITE(40,*) "Erro relativo =",erro_relativo*100,"%"
            WRITE(40,*) "Momento Fletor X máximo por Superposição de ",mxmax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
            WRITE(40,*) "Momento Fletor Y máximo por Superposição de ",mymax_calc/1000,"kN.m em X=",xmax_calc," Y=",ymax_calc
        CLOSE (40)
    END IF
END IF

IF(largura == comprimento) THEN
    k = 0
    DO i = 1, qtd_bordas
        IF(tipo_borda(i) == 'BEx' .OR. tipo_borda(i) == 'BEy') THEN
            k = k+1
        END IF
    END DO
    IF(k == 4) THEN
        WRITE(*,*) "teste:SOLUÇÃO DE RITZ"
        xmax_calc = 0.5*comprimento
        ymax_calc = 0.5*largura
        wmax_calc = 0.00126*carga*(comprimento**4)/mod_rigidez
        WRITE(*,*) "Deslocamento máximo por Navier de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
        erro_relativo = ABS(wmax_mdf-wmax_calc)/wmax_calc
        WRITE(*,*) "Erro relativo =",erro_relativo*100,"%"
        OPEN(UNIT=40, FILE='Deslocamento Maximo.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
            WRITE(40,*) "Deslocamento máximo por MDF de ",wmax_mdf*1000,"mm em X=",xmax_mdf," Y=",ymax_mdf
            WRITE(40,*) "Deslocamento máximo por Navier de ",wmax_calc*1000,"mm em X=",xmax_calc," Y=",ymax_calc
            WRITE(40,*) "Erro relativo =",erro_relativo*100,"%"
        CLOSE (40)
    END IF
END IF

tipo_grafico = 'Z'
DO WHILE (tipo_grafico /= 'S')
    WRITE(*,*) " "
    WRITE(*,*) "> Para gerar um gráfico, digite a opção:"
    WRITE(*,*) "D = Deslocamento;"
    WRITE(*,*) "R = Rotação;"
    WRITE(*,*) "C = Curvatura"
    WRITE(*,*) "CT = Curvatura Torsional"
    WRITE(*,*) "F = Momento Fletor"
    WRITE(*,*) "T = Momento Torsor"
    WRITE(*,*) "CO = Esforço Cortante"
    WRITE(*,*) "TE = Esforço Transverso Efetivo"
    WRITE(*,*) "FC = Força de Canto"
    WRITE(*,*) " "
    WRITE(*,*) "> Para voltar ao início, digite 'V'."
    WRITE(*,*) "> Para sair, digite 'S'."
    678 CONTINUE
    READ(*,*) tipo_grafico
    IF(tipo_grafico == 'S') THEN
        STOP
    ELSE IF(tipo_grafico == 'V') THEN
        DEALLOCATE(tipo_borda)
        DEALLOCATE(posicao_borda)
        DEALLOCATE(mat_pontos)
        DEALLOCATE(mat_coeficientes)
        DEALLOCATE(ipiv)
        DEALLOCATE(vet_cargas)
        DEALLOCATE(vetor_ordenada)
        DEALLOCATE(vetor_abscissa)
        DEALLOCATE(mat_deslocamento)
        DEALLOCATE(mat_expandida)
        DEALLOCATE(mat_fx)
        DEALLOCATE(mat_fy)
        GO TO 909
    ELSE IF(tipo_grafico /= 'D' .AND. tipo_grafico /= 'R' .AND. tipo_grafico /= 'C' .AND. tipo_grafico /= 'CT' .AND. &
            tipo_grafico /= 'F' .AND. tipo_grafico /= 'T' .AND. tipo_grafico /= 'CO' .AND. tipo_grafico /= 'TE' .AND. &
            tipo_grafico /= 'FC') THEN
        WRITE(*,*) "Opção inválida."
        GO TO 678
    END IF
    IF(tipo_grafico == 'R' .OR. tipo_grafico == 'C' .OR. tipo_grafico == 'F' .OR. tipo_grafico == 'CO' .OR. &
       tipo_grafico == 'TE') THEN
        WRITE(*,*) "Digite 'X' para a direção x ou 'Y' para a direção y."
        679 CONTINUE
        READ(*,*) direcao
        IF(direcao == 'S') THEN
            STOP
        ELSE IF(tipo_arquivo == 'V') THEN
            DEALLOCATE(tipo_borda)
            DEALLOCATE(posicao_borda)
            DEALLOCATE(mat_pontos)
            DEALLOCATE(mat_coeficientes)
            DEALLOCATE(ipiv)
            DEALLOCATE(vet_cargas)
            DEALLOCATE(vetor_ordenada)
            DEALLOCATE(vetor_abscissa)
            DEALLOCATE(mat_deslocamento)
            DEALLOCATE(mat_expandida)
            DEALLOCATE(mat_fx)
            DEALLOCATE(mat_fy)
            GO TO 909
        ELSE IF(direcao /= 'X' .AND. direcao /= 'Y') THEN
            WRITE(*,*) "Opção inválida."
            GO TO 679
        END IF
    END IF
    WRITE(*,*) "Digite 'I' para gerar um gráfico interativo ou 'P' para salvar em imagem .png."
    680 CONTINUE
    READ(*,*) tipo_arquivo
    IF(tipo_arquivo == 'S') THEN
        STOP
    ELSE IF(tipo_arquivo == 'V') THEN
        DEALLOCATE(tipo_borda)
        DEALLOCATE(posicao_borda)
        DEALLOCATE(mat_pontos)
        DEALLOCATE(mat_coeficientes)
        DEALLOCATE(ipiv)
        DEALLOCATE(vet_cargas)
        DEALLOCATE(vetor_ordenada)
        DEALLOCATE(vetor_abscissa)
        DEALLOCATE(mat_deslocamento)
        DEALLOCATE(mat_expandida)
        DEALLOCATE(mat_fx)
        DEALLOCATE(mat_fy)
        GO TO 909
    ELSE IF(tipo_arquivo /= 'I' .AND. tipo_arquivo /= 'P') THEN
        WRITE(*,*) "Opção inválida."
        GO TO 680
    END IF
    WRITE(*,*) "> Feche o terminal ou imagem para retornar às opções."
    IF(tipo_arquivo == 'I') THEN
        IF(tipo_grafico == 'D') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p plotD.plt')
        ELSE IF(tipo_grafico == 'CT') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p plotCT.plt')
        ELSE IF(tipo_grafico == 'T') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p plotT.plt')
        ELSE IF(tipo_grafico == 'FC') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p plotFC.plt')
        ELSE IF(direcao == 'X') THEN
            IF(tipo_grafico == 'R') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotRx.plt')
            ELSE IF(tipo_grafico == 'C') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotCx.plt')
            ELSE IF(tipo_grafico == 'F') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotFx.plt')
            ELSE IF(tipo_grafico == 'CO') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotCOx.plt')
            ELSE IF(tipo_grafico == 'TE') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotTEx.plt')
            END IF
        ELSE IF(direcao == 'Y') THEN
            IF(tipo_grafico == 'R') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotRy.plt')
            ELSE IF(tipo_grafico == 'C') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotCy.plt')
            ELSE IF(tipo_grafico == 'F') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotFy.plt')
            ELSE IF(tipo_grafico == 'CO') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotCOy.plt')
            ELSE IF(tipo_grafico == 'TE') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p plotTEy.plt')
            END IF
        END IF
    ELSE IF(tipo_arquivo == 'P') THEN
        IF(tipo_grafico == 'D') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p pngD.plt')
            CALL EXECUTE_COMMAND_LINE('graf-desloc.png')
        ELSE IF(tipo_grafico == 'CT') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p pngCT.plt')
            CALL EXECUTE_COMMAND_LINE('graf-curv-tors.png')
        ELSE IF(tipo_grafico == 'T') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p pngT.plt')
            CALL EXECUTE_COMMAND_LINE('graf-mom-tors.png')
        ELSE IF(tipo_grafico == 'FC') THEN
            CALL EXECUTE_COMMAND_LINE('gnuplot -p pngFC.plt')
            CALL EXECUTE_COMMAND_LINE('graf-for-canto.png')
        ELSE IF(direcao == 'X') THEN
            IF(tipo_grafico == 'R') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngRx.plt')
                CALL EXECUTE_COMMAND_LINE('graf-angulo-x.png')
            ELSE IF(tipo_grafico == 'C') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngCx.plt')
                CALL EXECUTE_COMMAND_LINE('graf-curv-x.png')
            ELSE IF(tipo_grafico == 'F') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngFx.plt')
                CALL EXECUTE_COMMAND_LINE('graf-mom-flet-x.png')
            ELSE IF(tipo_grafico == 'CO') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngCOx.plt')
                CALL EXECUTE_COMMAND_LINE('graf-cort-x.png')
            ELSE IF(tipo_grafico == 'TE') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngTEx.plt')
                CALL EXECUTE_COMMAND_LINE('graf-transv-x.png')
            END IF
        ELSE IF(direcao == 'Y') THEN
            IF(tipo_grafico == 'R') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngRy.plt')
                CALL EXECUTE_COMMAND_LINE('graf-angulo-y.png')
            ELSE IF(tipo_grafico == 'C') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngCy.plt')
                CALL EXECUTE_COMMAND_LINE('graf-curv-y.png')
            ELSE IF(tipo_grafico == 'F') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngFy.plt')
                CALL EXECUTE_COMMAND_LINE('graf-mom-flet-y.png')
            ELSE IF(tipo_grafico == 'CO') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngCOy.plt')
                CALL EXECUTE_COMMAND_LINE('graf-cort-y.png')
            ELSE IF(tipo_grafico == 'TE') THEN
                CALL EXECUTE_COMMAND_LINE('gnuplot -p pngTEy.plt')
                CALL EXECUTE_COMMAND_LINE('graf-transv-y.png')
            END IF
        END IF
    END IF
END DO

CONTAINS !Subrotinas
! Cada linha i tem todos os coeficientes dos pontos da placa para a equação aplicada no ponto k
! Cada coluna possui os coeficientes de um ponto presentes em todas as equações
SUBROUTINE GERMAIN_LAGRANGE
    ! WRITE(*,*) "teste:APLICAÇÃO DA EQUAÇÃO DE PLACAS DELGADAS SUJEITAS A PEQUENAS DEFLEXÕES"
    pt_i = j-2*x_div-8 !ponto i
    pt_jme1 = j-x_div-5 !ponto j-1
    pt_j = j-x_div-4 !ponto j
    pt_jma1 = j-x_div-3 !ponto j+1
    pt_kme2 = j-2 !ponto k-2
    pt_kme1 = j-1 !ponto k-1
    pt_k = j !ponto k
    pt_kma1 = j+1 !ponto k+1
    pt_kma2 = j+2 !ponto k+2
    pt_lme1 = j+x_div+3 !ponto l-1
    pt_l = j+x_div+4 !ponto l
    pt_lma1 = j+x_div+5 !ponto l+1
    pt_m = j+2*x_div+8 !ponto m
    mat_pontos(i,pt_i) = 1/(delta_y**4) !ponto i
    mat_pontos(i,pt_jme1) = 2/((delta_x*delta_y)**2) !ponto j-1
    mat_pontos(i,pt_j) = -4/(delta_y**4) - 4/((delta_x*delta_y)**2) !ponto j
    mat_pontos(i,pt_jma1) = 2/((delta_x*delta_y)**2) !ponto j+1
    mat_pontos(i,pt_kme2) = 1/(delta_x**4) !ponto k-2
    mat_pontos(i,pt_kme1) = -4/(delta_x**4) - 4/((delta_x*delta_y)**2)!ponto k-1
    mat_pontos(i,pt_k) = 6/(delta_x**4) + 6/(delta_y**4) + 8/((delta_x*delta_y)**2) !ponto k
    mat_pontos(i,pt_kma1) = -4/(delta_x**4) - 4/((delta_x*delta_y)**2) !ponto k+1
    mat_pontos(i,pt_kma2) = 1/(delta_x**4) !ponto k+2
    mat_pontos(i,pt_lme1) = 2/((delta_x*delta_y)**2) !ponto l-1
    mat_pontos(i,pt_l) = -4/(delta_y**4) - 4/((delta_x*delta_y)**2) !ponto l
    mat_pontos(i,pt_lma1) = 2/((delta_x*delta_y)**2) !ponto l+1
    mat_pontos(i,pt_m) = 1/(delta_y**4) !ponto m
END SUBROUTINE

SUBROUTINE BSA
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DA BORDA SIMPLESMENTE APOIADA"
    IF(posicao_borda(k) == 'SUP') THEN !Borda superior
        j = 3*x_div+15
        DO i = x_div+1, 2*x_div !Condição My=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  - mat_pontos(i,pt_i)*poisson*((delta_y/delta_x)**2)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + mat_pontos(i,pt_i)*(2+2*poisson*((delta_y/delta_x)**2))
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  - mat_pontos(i,pt_i)*poisson*((delta_y/delta_x)**2)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_i)
            mat_pontos(i,pt_i) = 0 !ponto i
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Borda direita
        j = 3*x_div+9
        DO i = x_div-1, pontos-1, x_div !Condição Mx=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  - mat_pontos(i,pt_kma2)*poisson*((delta_x/delta_y)**2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + mat_pontos(i,pt_kma2)*(2+2*poisson*((delta_x/delta_y)**2))
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  - mat_pontos(i,pt_kma2)*poisson*((delta_x/delta_y)**2)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma2) = 0 !ponto k+2
            j = j+x_div+4
        END DO
    ELSE IF(posicao_borda(k) == 'INF') THEN !Borda inferior
        j = pontos_extra-4*x_div-13
        DO i = pontos-2*x_div+1, pontos-x_div !Condição My=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  - mat_pontos(i,pt_m)*poisson*((delta_y/delta_x)**2)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + mat_pontos(i,pt_m)*(2+2*poisson*((delta_y/delta_x)**2))
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  - mat_pontos(i,pt_m)*poisson*((delta_y/delta_x)**2)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_m)
            mat_pontos(i,pt_m) = 0 !ponto m
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Borda esquerda
        j = 2*x_div+12
        DO i = 2, pontos-x_div+2, x_div !Condição Mx=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  - mat_pontos(i,pt_kme2)*poisson*((delta_x/delta_y)**2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + mat_pontos(i,pt_kme2)*(2+2*poisson*((delta_x/delta_y)**2))
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  - mat_pontos(i,pt_kme2)*poisson*((delta_x/delta_y)**2)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme2) = 0 !ponto k-2
            j = j+x_div+4
        END DO
    END IF
END SUBROUTINE

SUBROUTINE BE
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DA BORDA ENGASTADA"
    IF(posicao_borda(k) == 'SUP') THEN !Borda superior
        j = 3*x_div+15
        DO i = x_div+1, 2*x_div !Condição thetaY=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               + mat_pontos(i,pt_i)
            mat_pontos(i,pt_i) = 0 !ponto i
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Borda direita
        j = 3*x_div+9
        DO i = x_div-1, pontos-1, x_div !Condição thetaX=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               + mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma2) = 0 !ponto k+2
            j = j+x_div+4
        END DO
    ELSE IF(posicao_borda(k) == 'INF') THEN !Borda inferior
        j = pontos_extra-4*x_div-13
        DO i = pontos-2*x_div+1, pontos-x_div !Condição thetaY=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               + mat_pontos(i,pt_m)
            mat_pontos(i,pt_m) = 0 !ponto m
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Borda esquerda
        j = 2*x_div+12
        DO i = 2, pontos-x_div+2, x_div !Condição thetaX=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               + mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme2) = 0 !ponto k-2
            j = j+x_div+4
        END DO
    END IF
END SUBROUTINE


SUBROUTINE BD
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DA BORDA DESLIZANTE"
    IF(posicao_borda(k) == 'SUP') THEN !Borda superior
        j = 2*x_div+11
        DO i = 1, x_div !Condição Vy=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  - (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  - (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + (2+(4-2*poisson)*((delta_y/delta_x)**2))*mat_pontos(i,pt_i)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + (-2+(2*poisson-4)*((delta_y/delta_x)**2))*mat_pontos(i,pt_i)
            mat_pontos(i,pt_m) = mat_pontos(i,pt_m) & !ponto m
                               + mat_pontos(i,pt_i)
            mat_pontos(i,pt_i) = 0 !ponto i
            IF(i /= 1) THEN !Condição thetaY=0
                mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                      + mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_jme1) = 0 !ponto j-1
            END IF
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + mat_pontos(i,pt_j)
            mat_pontos(i,pt_j) = 0 !ponto j
            IF(i /= x_div) THEN
                mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                      + mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_jma1) = 0 !ponto j+1
            END IF
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Borda direita
        j = 3*x_div+10
        DO i = x_div, pontos, x_div !Condição Vx=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (2-poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  + (2-poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kme2) = mat_pontos(i,pt_kme2) & !ponto k-2
                                 + mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                              + (-2+(2*poisson-4)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + (2+(4-2*poisson)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma2) = 0 !ponto k+2
            IF(i /= x_div) THEN !Condição thetaX=0
                mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                      + mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_jma1) = 0 !ponto j+1
            END IF
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + mat_pontos(i,pt_kma1)
            mat_pontos(i,pt_kma1) = 0 !ponto k+1
            IF(i /= pontos) THEN
                mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                      + mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_lma1) = 0 !ponto l+1
            END IF
            j = j+x_div+4
        END DO
    ELSE IF(posicao_borda(k) == 'INF') THEN !Borda inferior
        j = pontos_extra-3*x_div-9
        DO i = pontos-x_div+1, pontos !Condição Vy=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  - (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  - (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + (-2+(2*poisson-4)*((delta_y/delta_x)**2))*mat_pontos(i,pt_m)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + (2+(4-2*poisson)*((delta_y/delta_x)**2))*mat_pontos(i,pt_m)
            mat_pontos(i,pt_i) = mat_pontos(i,pt_i) & !ponto i
                               + mat_pontos(i,pt_m)
            mat_pontos(i,pt_m) = 0 !ponto m
            !Condição thetaY=0
            IF(i /= pontos-x_div+1) THEN
                mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                      + mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_lme1) = 0 !ponto l-1
            END IF
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + mat_pontos(i,pt_l)
            mat_pontos(i,pt_l) = 0 !ponto l
            IF(i /= pontos) THEN
                mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                      + mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_lma1) = 0 !ponto l+1
            END IF
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Borda esquerda
        j = 2*x_div+11
        DO i = 1, pontos-x_div+1, x_div !Condição Vx=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (2+poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  + (2+poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kma2) = mat_pontos(i,pt_kma2) & !ponto k+2
                                  + mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + (2+(4-2*poisson)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + (-2+(2*poisson-4)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme2) = 0 !ponto k-2
            !Condição thetaX=0
            IF(i /= 1) THEN
                mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                      + mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_jme1) = 0 !ponto j-1
            END IF
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + mat_pontos(i,pt_kme1)
            mat_pontos(i,pt_kme1) = 0 !ponto k-1
            IF(i /= pontos-x_div+1) THEN
                mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                      + mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_lme1) = 0 !ponto l-1
            END IF
            j = j+x_div+4
        END DO
    END IF
END SUBROUTINE

SUBROUTINE BL
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DA BORDA LIVRE"
    IF(posicao_borda(k) == 'SUP') THEN !Borda superior
        j = 2*x_div+11
        DO i = 1, x_div !Condição Vy=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  - (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  - (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_i)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + (2+(4-2*poisson)*((delta_y/delta_x)**2))*mat_pontos(i,pt_i)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + (-2+(2*poisson-4)*((delta_y/delta_x)**2))*mat_pontos(i,pt_i)
            mat_pontos(i,pt_m) = mat_pontos(i,pt_m) & !ponto m
                               + mat_pontos(i,pt_i)
            mat_pontos(i,pt_i) = 0 !ponto i
            !Condição My=0
            IF(i /= 1) THEN
                mat_pontos(i,pt_kme2) = mat_pontos(i,pt_kme2) & !ponto k-2
                                      - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                      + (2+2*poisson*((delta_y/delta_x)**2))*mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                      - mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_jme1) = 0 !ponto j-1
            END IF
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_j)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               + (2+2*poisson*((delta_y/delta_x)**2))*mat_pontos(i,pt_j)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_j)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               - mat_pontos(i,pt_j)
            mat_pontos(i,pt_j) = 0 !ponto j
            IF(i /= x_div) THEN
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                      + (2+2*poisson*((delta_y/delta_x)**2))*mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_kma2) = mat_pontos(i,pt_kma2) & !ponto k+2
                                      - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                      - mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_jma1) = 0 !ponto j+1
            END IF
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Borda direita
        j = 3*x_div+10
        DO i = x_div, pontos, x_div !Condição theta=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (2-poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  + (2-poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kme2) = mat_pontos(i,pt_kme2) & !ponto k-2
                                  + mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + (-2+(2*poisson-4)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + (2+(4-2*poisson)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma2) = 0 !ponto k+2
            !Condição Mx=0
            IF(i /= x_div) THEN
                mat_pontos(i,pt_i) = mat_pontos(i,pt_i) & !ponto i
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                                   +(2+2*poisson*((delta_x/delta_y)**2))*mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                      - mat_pontos(i,pt_jma1)
                mat_pontos(i,pt_jma1) = 0 !ponto j+1
            END IF
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma1)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               +(2+2*poisson*((delta_x/delta_y)**2))*mat_pontos(i,pt_kma1)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_kma1)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  - mat_pontos(i,pt_kma1)
            mat_pontos(i,pt_kma1) = 0 !ponto k+1
            IF(i /= pontos) THEN
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                                   +(2+2*poisson*((delta_x/delta_y)**2))*mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_m) = mat_pontos(i,pt_m) & !ponto m
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                      - mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_lma1) = 0 !ponto l+1
            END IF
            j = j+x_div+4
        END DO
    ELSE IF(posicao_borda(k) == 'INF') THEN !Borda inferior
        j = pontos_extra-3*x_div-9
        DO i = pontos-x_div+1, pontos !Condição theta=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (poisson-2)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  - (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  - (2-poisson)*((delta_y/delta_x)**2)*mat_pontos(i,pt_m)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + (-2+(2*poisson-4)*((delta_y/delta_x)**2))*mat_pontos(i,pt_m)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + (2+(4-2*poisson)*((delta_y/delta_x)**2))*mat_pontos(i,pt_m)
            mat_pontos(i,pt_i) = mat_pontos(i,pt_i) & !ponto i
                               + mat_pontos(i,pt_m)
            mat_pontos(i,pt_m) = 0 !ponto m
            !Condição My=0
            IF(i /= pontos-x_div+1) THEN
                mat_pontos(i,pt_kme2) = mat_pontos(i,pt_kme2) & !ponto k-2
                                      - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                      + (2+2*poisson*((delta_y/delta_x)**2))*mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                      - mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_lme1) = 0 !ponto l-1
            END IF
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_l)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               + (2+2*poisson*((delta_y/delta_x)**2))*mat_pontos(i,pt_l)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_l)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               - mat_pontos(i,pt_l)
            mat_pontos(i,pt_l) = 0 !ponto l
            IF(i /= pontos) THEN
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                      + (2+2*poisson*((delta_y/delta_x)**2))*mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_kma2) = mat_pontos(i,pt_kma2) & !ponto k+2
                                      - poisson*((delta_y/delta_x)**2)*mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                      - mat_pontos(i,pt_lma1)
                mat_pontos(i,pt_lma1) = 0 !ponto l+1
            END IF
            j = j+1
        END DO
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Borda esquerda
        j = 2*x_div+11
        DO i = 1, pontos-x_div+1, x_div !Condição theta=0
            pt_i = j-2*x_div-8 !ponto i
            pt_jme1 = j-x_div-5 !ponto j-1
            pt_j = j-x_div-4 !ponto j
            pt_jma1 = j-x_div-3 !ponto j+1
            pt_kme2 = j-2 !ponto k-2
            pt_kme1 = j-1 !ponto k-1
            pt_k = j !ponto k
            pt_kma1 = j+1 !ponto k+1
            pt_kma2 = j+2 !ponto k+2
            pt_lme1 = j+x_div+3 !ponto l-1
            pt_l = j+x_div+4 !ponto l
            pt_lma1 = j+x_div+5 !ponto l+1
            pt_m = j+2*x_div+8 !ponto m
            mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                  + (2+poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                  + (2+poisson)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                                  + (poisson-2)*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kma2) = mat_pontos(i,pt_kma2) & !ponto k+2
                                  + mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + (2+(4-2*poisson)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + (-2+(2*poisson-4)*((delta_x/delta_y)**2))*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme2) = 0 !ponto k-2
            !Condição Mx=0
            IF(i/= 1) THEN
                mat_pontos(i,pt_i) = mat_pontos(i,pt_i) & !ponto i
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                                   +(2+2*poisson*((delta_x/delta_y)**2))*mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                                      - mat_pontos(i,pt_jme1)
                mat_pontos(i,pt_jme1) = 0 !ponto j-1
            END IF
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme1)
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               +(2+2*poisson*((delta_x/delta_y)**2))*mat_pontos(i,pt_kme1)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_kme1)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  - mat_pontos(i,pt_kme1)
            mat_pontos(i,pt_kme1) = 0 !ponto k-1
            IF(i /= pontos-x_div+1) THEN
                mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                                   +(2+2*poisson*((delta_x/delta_y)**2))*mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_m) = mat_pontos(i,pt_m) & !ponto m
                                   - poisson*((delta_x/delta_y)**2)*mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                                      - mat_pontos(i,pt_lme1)
                mat_pontos(i,pt_lme1) = 0 !ponto l-1
            END IF
            j = j+x_div+4
        END DO
    END IF
END SUBROUTINE

SUBROUTINE CANTO_BSABSA_BEBE_BSABE
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DO CANTO BSA-BSA, BE-BE e BSA-BE"
    IF(posicao_borda(k) == 'SUP') THEN !Canto superior direito
        DO i = 1, pontos !Condição w=0 superior
            DO j = 2*x_div+11, 3*x_div+10
                mat_pontos(i,j) = 0
            END DO
        END DO
        DO i = 1, x_div !Condição w=0 superior
            DO j = 1, pontos_extra
                mat_pontos(i,j) = 0
            END DO
        END DO
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Canto inferior direito
        DO i = 1, pontos !Condição w=0 direita
            DO j = 3*x_div+10, pontos_extra-2*x_div-10, x_div+4
                mat_pontos(i,j) = 0
            END DO
        END DO
        DO i = x_div, pontos, x_div !Condição w=0 direita
            DO j = 1, pontos_extra
                mat_pontos(i,j) = 0
            END DO
        END DO
    ELSE IF(posicao_borda(k) == 'INF') THEN !Canto inferior esquerdo
        DO i = 1, pontos !Condição w=0 inferior
            DO j = pontos_extra-3*x_div-10, pontos_extra-2*x_div-10
                mat_pontos(i,j) = 0
            END DO
        END DO
        DO i = pontos-x_div+1, pontos !Condição w=0 inferior
            DO j = 1, pontos_extra
                mat_pontos(i,j) = 0
            END DO
        END DO
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Canto superior esquerdo
        DO i = 1, pontos !Condição w=0 esquerda
            DO j = 2*x_div+11, pontos_extra-3*x_div-9, x_div+4
                mat_pontos(i,j) = 0
            END DO
        END DO
        DO i = 1, pontos-x_div+1, x_div !Condição w=0 esquerda
            DO j = 1, pontos_extra
                mat_pontos(i,j) = 0
            END DO
        END DO
    END IF
END SUBROUTINE

SUBROUTINE CANTO_BSABD
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DO CANTO BSA-BD"
    IF(posicao_borda(k) == 'SUP') THEN !Canto superior direito
        IF(tipo_borda(k) == 'BDx') THEN !Condição thetaY=0
            CALL BD
        ELSE IF(tipo_borda(k) == 'BSAx') THEN
            DO i = 1, pontos !Condição w=0 superior
                DO j = 2*x_div+11, 3*x_div+10
                    mat_pontos(i,j) = 0
                END DO
            END DO
            DO i = 1, x_div !Condição w=0 superior
                DO j = 1, pontos_extra
                    mat_pontos(i,j) = 0
                END DO
            END DO
        END IF
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Canto inferior direito
        IF(tipo_borda(k) == 'BDy') THEN !Condição thetaX=0
            CALL BD
        ELSE IF(tipo_borda(k) == 'BSAy') THEN
            DO i = 1, pontos !Condição w=0 direita
                DO j = 3*x_div+10, pontos_extra-2*x_div-10, x_div+4
                    mat_pontos(i,j) = 0
                END DO
            END DO
            DO i = x_div, pontos, x_div !Condição w=0 direita
                DO j = 1, pontos_extra
                    mat_pontos(i,j) = 0
                END DO
            END DO
        END IF
    ELSE IF(posicao_borda(k) == 'INF') THEN !Canto inferior esquerdo
        IF(tipo_borda(k) == 'BDx') THEN
            CALL BD
        ELSE IF(tipo_borda(k) == 'BSAx') THEN
            DO i = 1, pontos !Condição w=0 inferior
                DO j = pontos_extra-3*x_div-10, pontos_extra-2*x_div-10
                    mat_pontos(i,j) = 0
                END DO
            END DO
            DO i = pontos-x_div+1, pontos !Condição w=0 inferior
                DO j = 1, pontos_extra
                    mat_pontos(i,j) = 0
                END DO
            END DO
        END IF
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Canto superior esquerdo
        IF(tipo_borda(k) == 'BSAy') THEN
                CALL BD
            DO i = 1, pontos !Condição w=0 esquerda
                DO j = 2*x_div+11, pontos_extra-3*x_div-9, x_div+4
                    mat_pontos(i,j) = 0
                END DO
            END DO
            DO i = 1, pontos-x_div+1, x_div !Condição w=0 esquerda
                DO j = 1, pontos_extra
                    mat_pontos(i,j) = 0
                END DO
            END DO
        END IF
    END IF
END SUBROUTINE

SUBROUTINE CANTO_BSABL
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DO CANTO BSA-BL"
    IF(posicao_borda(k) == 'SUP') THEN !Canto superior direito
        !Condição Fc=0
        i = x_div
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              - mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_jma1) = 0 !ponto j+1
        IF(tipo_borda(k) == 'BLx') THEN
            CALL BL !Condição M=0
            i = x_div-1
            j = 3*x_div+9
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + 2*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma2) = 0 !ponto k+2
        ELSE IF(tipo_borda(k) == 'BSAx') THEN
            k = k+1
            CALL BL !Condição M=0
            k = k-1
            i = 2*x_div
            j = 4*x_div+14
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + 2*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_i) = 0 !ponto i
        END IF
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Canto inferior direito
        !Condição Fc=0
        i = pontos
        j = pontos_extra-2*x_div-10
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              - mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lma1) = 0 !ponto l+1
        IF(tipo_borda(k) == 'BLy') THEN
            CALL BL !Condição M=0
            i = pontos-1
            j = pontos_extra-2*x_div-10
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma1) = mat_pontos(i,pt_kma1) & !ponto k+1
                                  + 2*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_kma2) = 0 !ponto k+2
        ELSE IF(tipo_borda(k) == 'BSAy') THEN
            k = k+1
            CALL BL !Condição M=0
            k = k-1
            i = pontos-x_div
            j = pontos_extra-3*x_div-14
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_m)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + 2*mat_pontos(i,pt_m)
            mat_pontos(i,pt_m) = 0 !ponto m
        END IF
    ELSE IF(posicao_borda(k) == 'INF') THEN !Canto inferior esquerdo
        i = pontos-x_div+1
        j = pontos_extra-3*x_div-9
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lme1) = 0 !ponto l-1
        IF(tipo_borda(k) == 'BLx') THEN
            CALL BL !Condição M=0
            i = pontos-x_div+2
            j = pontos_extra-3*x_div-8
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + 2*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme2) = 0 !ponto k-2
        ELSE IF(tipo_borda(k) == 'BSAx') THEN
            k = k+1
            CALL BL !Condição M=0
            k = k-1
            i = pontos-2*x_div+1
            j = pontos_extra-4*x_div-15
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_m)
            mat_pontos(i,pt_l) = mat_pontos(i,pt_l) & !ponto l
                               + 2*mat_pontos(i,pt_m)
            mat_pontos(i,pt_m) = 0 !ponto m
        END IF
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Canto superior esquerdo
        i = 1
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jme1) = 0 !ponto j-1
        IF(tipo_borda(k) == 'BLy') THEN
            CALL BL !Condição M=0
            i = x_div+1
            j = 3*x_div+15
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_j) = mat_pontos(i,pt_j) & !ponto j
                               + 2*mat_pontos(i,pt_kma2)
            mat_pontos(i,pt_i) = 0 !ponto i
        ELSE IF(tipo_borda(k) == 'BSAy') THEN
            k = 1
            CALL BL !Condição M=0
            k = 4
            i = 2
            j = 2*x_div+12
            mat_pontos(i,pt_k) = mat_pontos(i,pt_k) & !ponto k
                               - mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme1) = mat_pontos(i,pt_kme1) & !ponto k-1
                                  + 2*mat_pontos(i,pt_kme2)
            mat_pontos(i,pt_kme2) = 0 !ponto k-2
        END IF
    END IF
END SUBROUTINE

SUBROUTINE CANTO_BDBL
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DO CANTO BD-BL"
    IF(posicao_borda(k) == 'SUP') THEN !Canto superior direito
        !Condição Fc=0
        i = x_div
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              - mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_jma1) = 0 !ponto j+1
        IF(tipo_borda(k) == 'BLx') THEN
            CALL BL !Condição M=0
            k = k+1
            CALL BD !Condição theta=0
            k = k-1
        ELSE IF(tipo_borda(k) == 'BDx') THEN
            CALL BD !Condição theta=0
            k = k+1
            CALL BL !Condição M=0
            k = k-1
        END IF
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Canto inferior direito
        !Condição Fc=0
        i = pontos
        j = pontos_extra-2*x_div-10
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              - mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lma1) = 0 !ponto l+1
        IF(tipo_borda(k) == 'BLy') THEN
            CALL BL !Condição M=0
            k = k+1
            CALL BD !Condição theta=0
            k = k-1
        ELSE IF(tipo_borda(k) == 'BDy') THEN
            CALL BD !Condição theta=0
            k = k+1
            CALL BL !Condição M=0
            k = k-1
        END IF
    ELSE IF(posicao_borda(k) == 'INF') THEN !Canto inferior esquerdo
        i = pontos-x_div+1
        j = pontos_extra-3*x_div-9
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lme1) = 0 !ponto l-1
        IF(tipo_borda(k) == 'BLx') THEN
            CALL BL !Condição M=0
            k = k+1
            CALL BD !Condição theta=0
            k = k-1
        ELSE IF(tipo_borda(k) == 'BDx') THEN
            CALL BD !Condição theta=0
            k = k+1
            CALL BL !Condição M=0
            k = k-1
        END IF
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Canto superior esquerdo
        i = 1
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jme1) = 0 !ponto j-1
        IF(tipo_borda(k) == 'BLy') THEN
            CALL BL !Condição M=0
            k = k+1
            CALL BD !Condição theta=0
            k = k-1
        ELSE IF(tipo_borda(k) == 'BDy') THEN
            CALL BD !Condição theta=0
            k = k+1
            CALL BL !Condição M=0
            k = k-1
        END IF
    END IF
END SUBROUTINE

SUBROUTINE CANTO_BDBD
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DO CANTO BD-BD"
    IF(posicao_borda(k) == 'SUP') THEN !Canto superior direito
        !Condição Fc=0
        i = x_div
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              - mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_jma1) = 0 !ponto j+1
        CALL BD !Condição theta=0
        k = k+1
        CALL BD !Condição theta=0
        k = k-1
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Canto inferior direito
        !Condição Fc=0
        i = pontos
        j = pontos_extra-2*x_div-10
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              - mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lma1) = 0 !ponto l+1
        CALL BD !Condição theta=0
        k = k+1
        CALL BD !Condição theta=0
        k = k-1
    ELSE IF(posicao_borda(k) == 'INF') THEN !Canto inferior esquerdo
        i = pontos-x_div+1
        j = pontos_extra-3*x_div-9
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lme1) = 0 !ponto l-1
        CALL BD !Condição theta=0
        k = k+1
        CALL BD !Condição theta=0
        k = k-1
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Canto superior esquerdo
        i = 1
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jme1) = 0 !ponto j-1
        CALL BD !Condição theta=0
        k = k+1
        CALL BD !Condição theta=0
        k = k-1
    END IF
END SUBROUTINE

SUBROUTINE CANTO_BLBL
    WRITE(*,*) "teste:CONDIÇÕES DE CONTORNO DO CANTO BL-BL"
    IF(posicao_borda(k) == 'SUP') THEN !Canto superior direito
        !Condição Fc=0
        i = x_div
        j = 2*x_div+10
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              - mat_pontos(i,pt_jma1)
        mat_pontos(i,pt_jma1) = 0 !ponto j+1
        CALL BL !Condição M=0
        k = k+1
        CALL BL !Condição M=0
        k = k-1
    ELSE IF(posicao_borda(k) == 'DIR') THEN !Canto inferior direito
        !Condição Fc=0
        i = pontos
        j = pontos_extra-2*x_div-10
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              - mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_lma1)
        mat_pontos(i,pt_lma1) = 0 !ponto l+1
        CALL BL !Condição M=0
        k = k+1
        CALL BL !Condição M=0
        k = k-1
    ELSE IF(posicao_borda(k) == 'INF') THEN !Canto inferior esquerdo
        i = pontos-x_div+1
        j = pontos_extra-3*x_div-9
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_jme1) = mat_pontos(i,pt_jme1) & !ponto j-1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_lme1)
        mat_pontos(i,pt_lme1) = 0 !ponto l-1
        CALL BL !Condição M=0
        k = k+1
        CALL BL !Condição M=0
        k = k-1
    ELSE IF(posicao_borda(k) == 'ESQ') THEN !Canto superior esquerdo
        i = 1
        j = 2*x_div+11
        pt_i = j-2*x_div-8 !ponto i
        pt_jme1 = j-x_div-5 !ponto j-1
        pt_j = j-x_div-4 !ponto j
        pt_jma1 = j-x_div-3 !ponto j+1
        pt_kme2 = j-2 !ponto k-2
        pt_kme1 = j-1 !ponto k-1
        pt_k = j !ponto k
        pt_kma1 = j+1 !ponto k+1
        pt_kma2 = j+2 !ponto k+2
        pt_lme1 = j+x_div+3 !ponto l-1
        pt_l = j+x_div+4 !ponto l
        pt_lma1 = j+x_div+5 !ponto l+1
        pt_m = j+2*x_div+8 !ponto m
        mat_pontos(i,pt_lma1) = mat_pontos(i,pt_lma1) & !ponto l+1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jma1) = mat_pontos(i,pt_jma1) & !ponto j+1
                              - mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_lme1) = mat_pontos(i,pt_lme1) & !ponto l-1
                              + mat_pontos(i,pt_jme1)
        mat_pontos(i,pt_jme1) = 0 !ponto j-1
        CALL BL !Condição M=0
        k = k+1
        CALL BL !Condição M=0
        k = k-1
    END IF
END SUBROUTINE

SUBROUTINE EXPANSAO
    IF(tipo_borda(k) == 'BEx' .OR. tipo_borda(k) == 'BDx') THEN
        IF(posicao_borda(k) == 'SUP') THEN
            DO j = 3, x_div+2
                mat_expandida(2,j) = mat_expandida(4,j)
            END DO
        ELSE IF(posicao_borda(k) == 'INF') THEN
            DO j = 3, x_div+2
                mat_expandida(y_div+1,j) = mat_expandida(y_div-1,j)
            END DO
        END IF
    ELSE IF(tipo_borda(k) == 'BEy' .OR. tipo_borda(k) == 'BDy') THEN
        IF(posicao_borda(k) == 'DIR') THEN
            DO i = 3, y_div+2
                mat_expandida(i,x_div+3) = mat_expandida(i,x_div+1)
            END DO
        ELSE IF(posicao_borda(k) == 'ESQ') THEN
            DO i = 3, x_div+2
                mat_expandida(i,2) = mat_expandida(i,4)
            END DO
        END IF
    ELSE IF(tipo_borda(k) == 'BSAx') THEN
        IF(posicao_borda(k) == 'SUP') THEN
            DO j = 3, x_div+2
                mat_expandida(2,j) = mat_expandida(4,j)
            END DO
        ELSE IF(posicao_borda(k) == 'INF') THEN
            DO j = 3, x_div+2
                mat_expandida(y_div+1,j) = mat_expandida(y_div-1,j)
            END DO
        END IF
    ELSE IF(tipo_borda(k) == 'BSAy') THEN
        IF(posicao_borda(k) == 'DIR') THEN
            DO i = 3, y_div+2
                mat_expandida(i,x_div+3) = mat_expandida(i,x_div+1)
            END DO
        ELSE IF(posicao_borda(k) == 'ESQ') THEN
            DO i = 3, x_div+2
                mat_expandida(i,2) = mat_expandida(i,4)
            END DO
        END IF
    END IF
END SUBROUTINE

! SUBROUTINE A_x
    ! mat_pontos(i,j-1) = -1/(2*delta_x) !ponto k-1
    ! mat_pontos(i,j+1) = 1/(2*delta_x) !ponto k+1
! END SUBROUTINE

! SUBROUTINE A_y
    ! mat_pontos(i,j-x_div-4) = -1/(2*delta_y) !ponto j
    ! mat_pontos(i,j+x_div+4) = 1/(2*delta_y) !ponto l
! END SUBROUTINE

! SUBROUTINE C_x
    ! mat_pontos(i,j-1) = 1/(delta_x**2) !ponto k-1
    ! mat_pontos(i,j) = -2/(delta_x**2) !ponto k
    ! mat_pontos(i,j+1) = 1/(delta_x**2) !ponto k+1
! END SUBROUTINE

! SUBROUTINE C_y
    ! mat_pontos(i,j-x_div-4) = 1/(delta_y**2) !ponto j
    ! mat_pontos(i,j) = -2/(delta_y**2) !ponto k
    ! mat_pontos(i,j+x_div+4) = 1/(delta_y**2) !ponto l
! END SUBROUTINE

! SUBROUTINE CT
    ! mat_pontos(i,j-x_div-5) = 1/(4*delta_x*delta_y) !ponto j-1
    ! mat_pontos(i,j-x_div-3) = -1/(4*delta_x*delta_y) !ponto j+1
    ! mat_pontos(i,j+x_div+3) = -1/(4*delta_x*delta_y) !ponto l-1
    ! mat_pontos(i,j+x_div+5) = 1/(4*delta_x*delta_y) !ponto l+1
! END SUBROUTINE

SUBROUTINE F_x
    DO i = 3, y_div+2
        DO j = 3, x_div+2
            mat_fx(i-2,j-2) = -mat_expandida(i-1,j)*mod_rigidez*poisson/(delta_y**2) & !ponto j
                            - mat_expandida(i,j-1)*mod_rigidez/(delta_x**2) & !ponto k-1
                            + mat_expandida(i,j)*mod_rigidez*2*(1/(delta_x**2)+poisson/(delta_y**2)) & !ponto k
                            - mat_expandida(i,j+1)*mod_rigidez/(delta_x**2) & !ponto k+1
                            - mat_expandida(i+1,j)*mod_rigidez*poisson/(delta_y**2) !ponto l
        END DO
    END DO
    DO k = 1, qtd_bordas
        IF(tipo_borda(k) == 'BSAy' .OR. tipo_borda(k) == 'BLy') THEN
            IF(posicao_borda(k) == 'DIR') THEN
                DO i = 1, y_div
                    mat_fx(i,1) = 0
                END DO
            ELSE IF(posicao_borda(k) == 'ESQ') THEN
                DO i = 1, y_div
                    mat_fx(i,x_div) = 0
                END DO
            END IF
        END IF
    END DO
END SUBROUTINE

SUBROUTINE F_y
    DO i = 3, y_div+2
        DO j = 3, x_div+2
            mat_fy(i-2,j-2) = -mat_expandida(i-1,j)*mod_rigidez/(delta_y**2) & !ponto j
                            - mat_expandida(i,j-1)*mod_rigidez*poisson/(delta_x**2) & !ponto k-1
                            + mat_expandida(i,j)*mod_rigidez*2*(1/(delta_y**2)+poisson/(delta_x**2)) & !ponto k
                            - mat_expandida(i,j+1)*mod_rigidez*poisson/(delta_x**2) & !ponto k+1
                            - mat_expandida(i+1,j)*mod_rigidez/(delta_y**2) !ponto l
        END DO
    END DO
    DO k = 1, qtd_bordas
        IF(tipo_borda(k) == 'BSAx' .OR. tipo_borda(k) == 'BLx') THEN
            IF(posicao_borda(k) == 'SUP') THEN
                DO j = 1, x_div
                    mat_fy(1,j) = 0
                END DO
            ELSE IF(posicao_borda(k) == 'INF') THEN
                DO j = 1, x_div
                    mat_fy(y_div,j) = 0
                END DO
            END IF
        END IF
    END DO
END SUBROUTINE

! SUBROUTINE T
    ! mat_pontos(i,j-x_div-5) = mod_rigidez*(poisson-1)/(4*delta_x*delta_y) !ponto j-1
    ! mat_pontos(i,j-x_div-3) = -mod_rigidez*(poisson-1)/(4*delta_x*delta_y) !ponto j+1
    ! mat_pontos(i,j+x_div+3) = -mod_rigidez*(poisson-1)/(4*delta_x*delta_y) !ponto l-1
    ! mat_pontos(i,j+x_div+5) = mod_rigidez*(poisson-1)/(4*delta_x*delta_y) !ponto l+1
! END SUBROUTINE

! SUBROUTINE CO_x
    ! mat_pontos(i,j-x_div-5) = mod_rigidez/(2*delta_x*(delta_y**2)) !ponto j-1
    ! mat_pontos(i,j-x_div-3) = -mod_rigidez/(2*delta_x*(delta_y**2)) !ponto j+1
    ! mat_pontos(i,j-2) = -mod_rigidez/(2*(delta_x**3)) !ponto k-2
    ! mat_pontos(i,j-1) = -mod_rigidez*(1/(delta_x**3)+1/(delta_x*(delta_y**2))) !ponto k-1
    ! mat_pontos(i,j+1) = mod_rigidez*(1/(delta_x**3)+1/(delta_x*(delta_y**2))) !ponto k+1
    ! mat_pontos(i,j+2) = -mod_rigidez/(2*(delta_x**3)) !ponto k+2
    ! mat_pontos(i,j+x_div+3) = mod_rigidez/(2*delta_x*(delta_y**2)) !ponto l-1
    ! mat_pontos(i,j+x_div+5) = -mod_rigidez/(2*delta_x*(delta_y**2)) !ponto l+1
! END SUBROUTINE

! SUBROUTINE CO_y
    ! mat_pontos(i,j-2*x_div-8) = mod_rigidez/(2*(delta_y**3)) !ponto i
    ! mat_pontos(i,j-x_div-5) = mod_rigidez/(2*(delta_x**2)*delta_y) !ponto j-1
    ! mat_pontos(i,j-x_div-4) = -mod_rigidez*(1/((delta_x**2)*delta_y)+1/(delta_y**3)) !ponto j
    ! mat_pontos(i,j-x_div-3) = mod_rigidez/(2*(delta_x**2)*delta_y) !ponto j+1
    ! mat_pontos(i,j+x_div+3) = -mod_rigidez/(2*(delta_x**2)*delta_y) !ponto l-1
    ! mat_pontos(i,j+x_div+4) = mod_rigidez*(1/((delta_x**2)*delta_y)+1/(delta_y**3)) !ponto l
    ! mat_pontos(i,j+x_div+5) = -mod_rigidez/(2*(delta_x**2)*delta_y) !ponto l+1
    ! mat_pontos(i,j+2*x_div+8) = -mod_rigidez/(2*(delta_y**3)) !ponto m
! END SUBROUTINE

! SUBROUTINE TE_x
    ! mat_pontos(i,j-x_div-5) = mod_rigidez*(2-poisson)/(2*delta_x*(delta_y**2)) !ponto j-1
    ! mat_pontos(i,j-x_div-3) = mod_rigidez*(poisson-2)/(2*delta_x*(delta_y**2)) !ponto j+1
    ! mat_pontos(i,j-2) = -mod_rigidez/(2*(delta_x**3)) !ponto k-2
    ! mat_pontos(i,j-1) = -mod_rigidez*(1/(delta_x**3)+(2-poisson)/(delta_x*(delta_y**2))) !ponto k-1
    ! mat_pontos(i,j+1) = mod_rigidez*(1/(delta_x**3)+(poisson-2)/(delta_x*(delta_y**2))) !ponto k+1
    ! mat_pontos(i,j+2) = -mod_rigidez/(2*(delta_x**3)) !ponto k+2
    ! mat_pontos(i,j+x_div+3) = mod_rigidez*(2-poisson)/(2*delta_x*(delta_y**2)) !ponto l-1
    ! mat_pontos(i,j+x_div+5) = mod_rigidez*(poisson-2)/(2*delta_x*(delta_y**2)) !ponto l+1
! END SUBROUTINE

! SUBROUTINE TE_y
    ! mat_pontos(i,j-2*x_div-8) = mod_rigidez/(2*(delta_y**3)) !ponto i
    ! mat_pontos(i,j-x_div-5) = mod_rigidez*(2-poisson)/(2*(delta_x**2)*delta_y) !ponto j-1
    ! mat_pontos(i,j-x_div-4) = -mod_rigidez*((2-poisson)/((delta_x**2)*delta_y)+1/(delta_y**3)) !ponto j
    ! mat_pontos(i,j-x_div-3) = mod_rigidez*(2-poisson)/(2*(delta_x**2)*delta_y) !ponto j+1
    ! mat_pontos(i,j+x_div+3) = mod_rigidez*(poisson-2)/(2*(delta_x**2)*delta_y) !ponto l-1
    ! mat_pontos(i,j+x_div+4) = mod_rigidez*((2-poisson)/((delta_x**2)*delta_y)+1/(delta_y**3)) !ponto l
    ! mat_pontos(i,j+x_div+5) = mod_rigidez*(poisson-2)/(2*(delta_x**2)*delta_y) !ponto l+1
    ! mat_pontos(i,j+2*x_div+8) = -mod_rigidez/(2*(delta_y**3)) !ponto m
! END SUBROUTINE

! SUBROUTINE FC
    ! mat_aux(i,j) = &
    ! mat_pontos(i,j-x_div-5) = mod_rigidez*(poisson-1)/(2*delta_x*delta_y) !ponto j-1
    ! mat_pontos(i,j-x_div-3) = -mod_rigidez*(poisson-1)/(2*delta_x*delta_y) !ponto j+1
    ! mat_pontos(i,j+x_div+3) = -mod_rigidez*(poisson-1)/(2*delta_x*delta_y) !ponto l-1
    ! mat_pontos(i,j+x_div+5) = mod_rigidez*(poisson-1)/(2*delta_x*delta_y) !ponto l+1
! END SUBROUTINE

END PROGRAM PLACAS_DELGADAS