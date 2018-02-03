!
! Last revision by nelsonvn on 20160523
!
!!!Program DirFluxo4

subroutine fluxdir(nlinpix, ncolpix, dirh, areah, outlet, &
   nlincel, ncolcel, mask, dirl, outrow, outcol, plin, pcol, limiteinc, limitecam)


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!    Programa para gerar direcoes de fluxo em um malha de BAIXA resolucao a partir das
! direcoes de fluxo de uma malha de ALTA resolucao.
!    IMPORTANTE: as grades de baixa e alta resolucao devem ser multiplas (ex: 10km e 100m)
!
! Algoritmo adaptado do descrito no artigo de Reed (2003) e descrito em Paz et al.(2006).
!
! Considerando: grade de BAIXA resolucao = celulas
!               grade de ALTA resolucao = pixels
!     
!
! (i) O algoritmo identifica em cada celula o pixel exutorio, definido como aquele que drena
! a maior area de drenagem dentre todos os pixels contidos na celula, desde que atenda
! a outro criterio: o fluxo principal do escoamento a montante dele, dentro da celula, deve
! apresentar um comprimento minimo. Se nao satisfaz este ultimo criterio, avalia se ele eh o 
! pixel que drena a maior parte da celula; se for, eh aceito como pixel exutorio, caso contrario
! procura o seguinte em termos de area de drenagem acumulada para ser testado e assim por diante.
! (ii) Para cada celula, percorre-se o caminho do fluxo desde o pixel exutorio ateh encontrar
! o pixel exutorio de uma celula vizinha. Entao eh checado o incremento de area de drenagem e,
! se for superior ao minimo estabelecido, define-se a direcao da celula analisada. Se nao, 
! continua o percurso..
! (iii) Existem diversas situacoes especiais: ver artigos ou o codigo...
!
!
! Arquivos de entrada: 
! * arquivo com direcoes de fluxo de ALTA resolucao
! * arquivo com areas de drenagem acumuladas de ALTA resolucao (em km2)
! * (opcional) mascara da regiao de mar de BAIXA resolucao
!
!
! IMPORTANTE: 

! 1. O parametro LIMITEINC define a area de drenagem incremental minima para estabelecimento
! das direcoes de fluxo. Para celulas de 10 km x 10 km e pixels de cerca de 200 m, o valor de
! 100 km2 para esse parametro foi o mais indicado.
!
! 2. O parametro LIMITECAM define o tamanho do caminho minimo a montante do pixel testado 
! para que ele seja aceito como pixel exutorio. Valor indicado para esse parametro, considerando
! pixels em torno de 200 m e celulas de 10 km, eh de 1 a 2 km.
!
! 3. Opcionalmente, pode ser lida imagem raster com mascara definindo regiao de mar (ou outra sem
! interesse), para que tal regiao nao seja tratada como as demais pelo algoritmo.
! OPCAO definida pelo usuario na tela!!
!
! -------------------------
! Adriano Rolim da Paz - julho/2005
! adrianorpaz@yahoo.com.br
! VERSAO JULHO/2006
! VERSAO PARA DISTRIBUICAO PELA INTERNET ABR2008
!


!>>>>>>>>>>>>>>>>>>>> DEFINICAO DAS VARIAVEIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

implicit none

! Dummy arguments - for R compatibility
integer, intent(in) :: nlinpix, ncolpix, nlincel, ncolcel, plin, pcol, limitecam
real, intent(in) :: limiteinc
!integer*2, intent(in) :: dirh(ncolpix,nlinpix)
!integer*2, intent(in) :: mask(ncolcel,nlincel)
integer, intent(in) :: dirh(ncolpix,nlinpix)
integer, intent(in) :: mask(ncolcel,nlincel)
real, intent(in) :: areah(ncolpix,nlinpix)
!integer*2, intent(out) :: outlet(ncolpix,nlinpix)
!integer*2, intent(out) :: dirl(ncolcel,nlincel)
integer, intent(out) :: outlet(ncolpix,nlinpix)
integer, intent(out) :: dirl(ncolcel,nlincel)
integer, intent(out) :: outrow(ncolcel,nlincel)
integer, intent(out) :: outcol(ncolcel,nlincel)

! Local variables
integer :: i, j



! Original variables
!!!integer :: nlinpix, ncolpix, respmar
integer :: linpix, colpix
!integer*2 :: dirpix(nlinpix, ncolpix)
!integer*2 :: exu(nlinpix, ncolpix)
integer :: dirpix(nlinpix, ncolpix)
integer :: exu(nlinpix, ncolpix)

!integer*2 :: dircel(nlincel, ncolcel)
!integer*2 :: marcel(nlincel, ncolcel)
!integer*2 :: vizmarcel(nlincel, ncolcel)
!integer*2 :: histocel(nlincel, ncolcel)
!integer*2 :: histo(nlinpix, ncolpix)
!integer*2 :: marcacel(nlincel, ncolcel)
integer :: dircel(nlincel, ncolcel)
integer :: marcel(nlincel, ncolcel)
integer :: vizmarcel(nlincel, ncolcel)
integer :: histocel(nlincel, ncolcel)
integer :: histo(nlinpix, ncolpix)
integer :: marcacel(nlincel, ncolcel)

integer :: linpixex(nlincel, ncolcel)
integer :: colpixex(nlincel, ncolcel)

real :: aacpix(nlinpix, ncolpix)
real :: areaaux(plin * pcol)



!!!character(len=5) :: texto1,texto2   
!!!character(len=9) :: texto3,texto4  
!!!character(len=8) :: texto5          
!!!character(len=12) :: texto6        
!!!integer :: pixnul   
!!!real :: coordxie,coordyie,respix,rescel     
!!!integer*2,allocatable:: dirpix(:,:),exu(:,:)   ! Flux direction and oulet grids
!!!real,allocatable:: aacpix(:,:)   ! Accumulated area grid
!!!character (len=40) :: texto7
!!!character (len=13) :: texto8

!!!integer :: ncolpix,nlinpix
!!!integer :: nlincel,ncolcel
!!!integer :: plin,pcol     
integer :: lincel, colcel, fim
!!!integer*2,allocatable:: dircel(:,:),marcel(:,:),vizmarcel(:,:),histocel(:,:),histo(:,:)
!!!integer*2,allocatable:: marcacel(:,:)
integer histomax

integer :: linpixini,colpixini,linpixfin,colpixfin
integer :: lincelini,colcelini,lincelfin,colcelfin

real MaiorArea
integer linexu,colexu
!!!integer,allocatable:: linpixex(:,:)
!!!integer,allocatable:: colpixex(:,:)

integer linaux,colaux
integer caminho
integer diraux
real auxl,auxc,auxl2,auxc2
integer lincelaux,colcelaux
!!!real AreaInc,LimiteInc
real AreaInc
integer deltalin,deltacol

integer foraviz,exut
!!!integer dd,i
integer dd
integer A,B,C,D,E,F,G,H
integer ArcView

!!!integer LimiteCam
!!!real LimiteCamkm,LimiteCamr
!!!real :: LimiteCamr

integer montante
integer foracel
real MaiorAreaExu
integer contacam
integer passo2
integer linaux2,colaux2
real MaiorArea2
integer linaux3,colaux3
!!!real,allocatable:: AreaAux(:)
integer ex,nex,auxex
real AreaJus,AreaMont

integer cruzamento
integer dircelSE,dircelSD,dircelIE,dircelID
integer correcao
integer aacpixSE,aacpixSD,aacpixIE,aacpixID

integer :: dlin(128),dcol(128),ddaux(8),inv(128)
integer laux,caux

integer dl, dc,invaux
integer lincelaux2,colcelaux2
integer caux99,laux99

integer cont,linexuH,colexuH

!!!real respalta,respbaixa
!!!real plinr,pcolr



! Copy INPUT data - for R compatibility
do j = 1, nlinpix
   do i = 1, ncolpix
      dirpix(j,i) = dirh(i,j)
   end do
end do
do j = 1, nlinpix
   do i = 1, ncolpix
      aacpix(j,i) = areah(i,j)
   end do
end do
do j = 1, nlincel
   do i = 1, ncolcel
      marcel(j,i) = mask(i,j)
   end do
end do



!! CABECALHO DO PROGRAMA 

!!!write(*,*)
!!!write(*,*)
!!!write(*,*) "-----------------------PROGRAMA DirFluxo----------------------------"
!!!write(*,*) "-                                                                  -"
!!!write(*,*) "- @ OBJETIVO: fazer upscaling de direcoes de fluxo                 -"
!!!write(*,*) "-                                                                  -"
!!!write(*,*) "- @@ ARQUIVOS DE ENTRADA:                                          -"
!!!write(*,*) "-  1.DirAlta.rst e .rdc (dir. de fluxo de alta resol./int/bin/IDR.)-"
!!!write(*,*) "-  2.AreaAlta.rst  (areas dren. acumuladas - raster/km2/real/bin/I)-"
!!!write(*,*) "-  3.MASCARA.rst    (int/bin/baixa resol., arq. raster do IDRISI)  -"
!!!write(*,*) "-  (obs: este eh opcional; pixels tem valor 0 ou 1, conforme       -"
!!!write(*,*) "-   representem regiao para entrar no calculo ou nao, respect)     -"
!!!write(*,*) "-                                                                  -"
!!!write(*,*) "- @@@ ARQUIVOS DE SAIDA:                                           -"
!!!write(*,*) "-  1.DirBaixa.rst  (dir. de fluxo de baixa resol./int/bin/IDRISI)  -"
!!!write(*,*) "-  2.PixelExu.dat  (matriz ascii com lin/col dos pixels exutorios) -"
!!!write(*,*) "-  3.PixelExu.rst  (raster alta res. 1(pixel exut.)/0 int/bin/IDR) -"
!!!write(*,*) "-                                                                  -" 
!!!write(*,*) "- @@@@ algoritmo descrito em Paz et al(2005) c/ base em Reed(2003) -"                                                     
!!!write(*,*) "-                                                                  -"
!!!write(*,*) "- @@@@@ contato: Adriano Rolim da Paz - adrianorpaz@yahoo.com.br   -"
!!!write(*,*) "-                Instituto de Pesquisas Hidraulicas (IPH/UFRGS)    -" 
!!!write(*,*) "-                                                                  -" 
!!!write(*,*) "- @@@@@@ versao: JULHO/2006                                        -" 
!!!write(*,*) "- @@@@@@ versao p/ distribuicao internet: ABR/2008                 -" 
!!!write(*,*) "--------------------------------------------------------------------"
!!!write(*,*) "(tecle enter)"
!!!read(*,*)


!!!!limpa a tela (apelando...)
!!!do i=1,200
!!!  write(*,*)
!!!end do



!>>>>>>>>>>>>>>>>>> DEFINICAO DOS PARAMETROS E LEITURA DOS ARQUIVOS DE ENTRADA >>>>>>>>>>>>>>>>
!
!!!write(*,*) "1. definindo parametros e lendo arquivos de entrada..."
!!!write(*,*)



!!!write(*,*) "Vai ler arquivo raster (inteiro/binario) de direcoes"
!!!write(*,*) "de fluxo de alta resolucao"
!!!write(*,*) "DirAlta.rdc e DirAlta.rst. CONTINUAR? (tecle enter)"
!!!read(*,*)
!!!write(*,*)



!ATENCAO PARA O VALOR NUMERICO DAS DIRECOES

!   G  H  A          ArcView:  32 64 128    MGB-IPH:  64  128  1 
!   F  *  B                    16  *  1               32   *   2
!   E  D  C                     8  4  2               16   8   4



!
! Read high resolution header (RDC) file
!
!!!   !leitura do arquivo .RDC das direcoes de fluxo de alta resolucao com numero de linhas e colunas
!!!	open(10,file='DIRALTA.rdc')
!!!	read(10,'(A)') texto7
!!!	read(10,'(A)') texto7
!!!	read(10,'(A)') texto7
!!!	read(10,'(A)') texto7
!!!	read(10,'(A,1I10)') texto8,ncolpix
!!!	read(10,'(A,1I10)') texto8,nlinpix
!!!	close(10)

    ! aloca variaveis em funcao do numero de linhas e colunas do arquivo .RDC
!!!  allocate(dirpix(nlinpix,ncolpix))
!!!	allocate(aacpix(nlinpix,ncolpix))
!!!	allocate(exu(nlinpix,ncolpix))

    ! leitura do arquivo raster .RST de direcoes de fluxo de alta resolucao - arquivo binario/inteiro
!!!	open(10,file='DIRALTA.rst',status='old',form='unformatted',access='direct',RECL=2*ncolpix)
!!!	do linpix=1,nlinpix
!!!	  read(10,REC=linpix) (dirpix(linpix,colpix),colpix=1,ncolpix)
!!!	end do
!!!	close(10)

!!!write(*,*) "ok"
!!!write(*,*)


!limpa a tela (apelando...)
!!!do i=1,200
!!!  write(*,*)
!!!end do

!!!write(*,*) "Vai ler arquivo raster (real/binario) de areas acumuladas de alta resolucao"
!!!write(*,*) "AreaAlta.rdc e AreaAlta.rst. CONTINUAR? (tecle enter)"
!!!read(*,*)

	! leitura do arquivo raster .RST de Areas Acumuladas de alta resolucao - arquivo binario/real
!!!	open(20,file='AREAALTA.rst',status='old',form='unformatted',access='direct',RECL=4*ncolpix)
!!!	do linpix=1,nlinpix
!!!	  read(20,REC=linpix) (aacpix(linpix,colpix),colpix=1,ncolpix)
!!!	end do
!!!	close(20)

!!!write(*,*) "ok"
!!!write(*,*)


!limpa a tela (apelando...)
!!!do i=1,200
!!!  write(*,*)
!!!end do

! definicao do numero de linhas e de colunas do arquivo com direcoes
! de fluxo de BAIXA resolucao
! (nlinpix e ncolpix devem ser multiplos deles)

!para ALTA resolucao de 0.1 grau e BAIXA resolucao de 0.005 grau (aprox 500m)
!-> reducao de 20x (cada celula contem 20 x 20 pixels)

!para ALTA resolucao de 0.1 grau e BAIXA resolucao de 0.0025 grau (aprox 250m)
!-> reducao de 40x (cada celula contem 40 x 40 pixels)

!!!write(*,*) "Informe resolucao da grade de ALTA resolucao (em graus)"
!!!read(*,*) respalta
!!!write(*,*) "Informe resolucao da grade de BAIXA resolucao (em graus)"
!!!read(*,*) respbaixa



!!!plinr=(respbaixa/respalta)
!!!plin=int(plinr)
!!!pcol=plin

!!!nlincel=nlinpix/plin
!!!ncolcel=ncolpix/pcol



!!!write(*,*) "Numero de linhas e colunas da grade de BAIXA resolucao: ",nlincel,ncolcel
!!!write(*,*) "continuar? (tecle enter). Se nao, cancele execucao..."
!!!read(*,*)
!!!write(*,*)

!!!allocate (dircel(nlincel,ncolcel))
!!!allocate (marcel(nlincel,ncolcel))
!!!allocate (vizmarcel(nlincel,ncolcel))
!!!allocate (linpixex(nlincel,ncolcel))
!!!allocate (colpixex(nlincel,ncolcel))
!!!allocate (areaaux(plin*pcol))

!!!allocate (histocel(nlincel,ncolcel))
!!!allocate (histo(nlinpix,ncolpix))
!!!allocate (marcacel(nlincel,ncolcel))

vizmarcel=0


!limpa a tela (apelando...)
!!!do i=1,200
!!!  write(*,*)
!!!end do

!!!write(*,*)
!!!write(*,*) "Na janela de trabalho existe regiao de mar (ou outra sem interesse)??"
!!!write(*,*) "nao -> nao precisa ler raster que define a regiao de mar"
!!!write(*,*) "sim -> precisa ler raster que define a regiao de mar (BAIXA resolucao)"
!!!write(*,*) "digite: 0 (nao vai leh raster); 1 (vai leh raster); 2 (parar execucao)"
!!!write(*,*)
!!!write(*,*) "nome do arquivo: MASCARA.RST. Arquivo binario/inteiro (BAIXA resolucao)"
!!!read(*,*) RespMar
!!!write(*,*)

!!!if (RespMar==1) then
	! leitura do arquivo contendo mascara que diferencia continente do mar - arquivo binario/inteiro
	!(mesmo numero de linhas e colunas (mesma resolucao) do arquivo de BAIXA resolucao do MNT, 
	! sendo valor 0 para regiao de mar e 1 para continente)
!!!	open(21,file='Mascara.rst',status='old',form='unformatted',access='direct',RECL=2*ncolcel)
!!!	do lincel=1,nlincel
!!!	  read(21,REC=lincel) (marcel(lincel,colcel),colcel=1,ncolcel)
!!!	end do
!!!	close(21)
!!!else
!!!  if (RespMar==0) then
!!!    marcel=0
!!!  else
!!!    stop
!!!  end if
!!!end if

!!!write(*,*) "ok"
!!!write(*,*)


!limpa a tela (apelando...)
!!!do i=1,200
!!!  write(*,*)
!!!end do

!limite minimo da AREA DRENAGEM ACUMULADA INCREMENTAL para definir pixel exutorio (EM KM2)!!!
!!!write(*,*) "Informe area de drenagem incremental minima (em km2)"
!!!write(*,*) "valor usual: 100 km2"
!!!read(*,*) LimiteInc
!!!write(*,*)

!limite minimo do CAMINHO DE MONTANTE para definir pixel exutorio (EM KM)!!!
!!!write(*,*) "Informe caminho minimo de montante (em km)"
!!!write(*,*) "valor usual entre 1 e 2 km" 
!!!read(*,*) LimiteCamkm
!!!LimiteCamR=LimiteCamkm/100.0/respalta
!!!LimiteCam=int(LimiteCamR)

!!!write(*,*) "ok"
!!!write(*,*)


!limpa a tela (apelando...)
!!!do i=1,200
!!!  write(*,*)
!!!end do



!definicao da numeracao das direcoes
A=1   
B=2   
C=4   
D=8   
E=16  
F=32  
G=64  
H=128 



!definicao do vetor de direcoes ddaux
ddaux(1)=1
ddaux(2)=2
ddaux(3)=4
ddaux(4)=8
ddaux(5)=16
ddaux(6)=32
ddaux(7)=64
ddaux(8)=128



!definicao da posicao relativa dos pixels vizinhos
dlin(A)=-1
dcol(A)=1
dlin(B)=0
dcol(B)=1
dlin(C)=1
dcol(C)=1
dlin(D)=1
dcol(D)=0
dlin(E)=1
dcol(E)=-1
dlin(F)=0
dcol(F)=-1
dlin(G)=-1
dcol(G)=-1
dlin(H)=-1
dcol(H)=0



!definicao das direcoes opostas
inv(A)=E
inv(B)=F
inv(C)=G
inv(D)=H
inv(E)=A
inv(F)=B
inv(G)=C
inv(H)=D



! inicializacao de alguns parametros
lincel=1
colcel=1
dircel=-99
linpixex=1
colpixex=1

foraviz=0

colpixini=1
colpixfin=pcol    
linpixini=1
linpixfin=plin

fim=0


histo=0
histocel=0

!!!write(*,*) "INICIO DO ALGORITMO...(aguarde conclusao)"
!!!write(*,*)

!>>>>>>>>>>>>>>>>>>>>>>>> DETERMINACAO DOS PIXELS EXUTORIOS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!!write(*,*) "2. procurando pixel exutorio em cada celula..."


! identificacao do pixel exutorio em cada celula
! (eh aquele cujo numero de pixels drenados eh superior a um valor estabelecido (LimitInc),
! satisfazendo tambem a condicao de que o numero de pixels a montante no caminho preferencial
! localizados dentro da celula seja maior do que um valor minimo (LimiteCam))

do while (fim==0)

  MaiorArea=0
  linexu=0
  colexu=0
  montante=0
  passo2=0
  exut=0
  AreaAux=0
	write(*,*) lincel,colcel
  if (marcel(lincel,colcel)==0) then
  do while (passo2==0)

  !loop para achar pixel de maior area de drenagem (procura apenas nas bordas da celula)
	! (1as e ultimas linhas e colunas)
	if (exut==1) then      !indica que algum pixel potencial exutorio encontrado foi descartado
	  nex=nex+1
	  AreaAux(nex)=MaiorAreaExu
      MaiorAreaExu=0
	else                   !indica que eh a primeira vez que procura um pixel exutorio
	  AreaAux=0
	  MaiorAreaExu=0
	  nex=0
	end if
		do linpix=linpixini,linpixfin
		  do colpix=colpixini,colpixfin
		   if ((linpix==linpixini).OR.(linpix==linpixfin).OR.(colpix==colpixini).OR.(colpix==colpixfin)) then
			   auxex=0
			   if (aacpix(linpix,colpix)>=MaiorAreaExu) then
				 if (nex>0) then
				   do ex=1,nex
					 if (aacpix(linpix,colpix)<AreaAux(ex)) then
					   auxex=auxex+1
					 end if
				   end do
				   if (auxex==nex) then
					   MaiorAreaExu=aacpix(linpix,colpix)
   					   linexu=linpix
					   colexu=colpix
					   exut=1
				   end if
				 else
				   MaiorAreaExu=aacpix(linpix,colpix)
   				   linexu=linpix
				   colexu=colpix
				   exut=1
				 end if	 
			   end if
			end if
		  end do
		end do
		!verificacao do criterio de numero minimo de pixels a montante no caminho do pixel exutorio
		linaux2=linexu
		colaux2=colexu
		contacam=0
		linaux=linexu
		colaux=colexu
		linaux3=linexu
		colaux3=colexu
		montante=0
		AreaJus=MaiorAreaExu
		do while (montante==0)
		  !determina pixel de montante e sua area drenagem  (eh o pixel vizinho de maior area, desde
		  !que area menor do que o do pixel jusante, que drena para ele)
		  MaiorArea=0
		  do linpix=linaux2-1,linaux2+1
			do colpix=colaux2-1,colaux2+1
			  if ((linpix>=1).AND.(linpix<=nlinpix).AND.(colpix>=1).AND.(colpix<=ncolpix)) then
				if ((linpix/=linaux2).OR.(colpix/=colaux2)) then
				  if ((linpix/=linaux3).OR.(colpix/=colaux3)) then
					if (aacpix(linpix,colpix)>MaiorArea) then
					  if (aacpix(linpix,colpix)<AreaJus) then
						!checa se ele drena para o pixel
						dl=linpix-linaux2
						dc=colpix-colaux2
						dd=dirpix(linpix,colpix)
						invaux=inv(dd)
						do i=1,8					 
						  if (invaux==ddaux(i)) then			   
							if ((dlin(invaux)==dl).AND.(dcol(invaux)==dc)) then					     
							  MaiorArea=aacpix(linpix,colpix)					
							  linaux=linpix
							  colaux=colpix								  
							end if
						  end if
						end do
					  end if
					end if
				  end if
				end if
			  end if
			end do
		  end do
		  AreaJus=MaiorArea
		  linaux3=linaux2
		  colaux3=colaux2
		  
		  !verifica se pixel de montante pertence a mesma celula
		  !    (linha)
		  if (linaux<=plin) then
  			lincelaux=1
		  else
			auxl=linaux/real(plin)
			auxl2=auxl-int(auxl)
			if (auxl2>0) then
			  lincelaux=int(auxl)+1
			else
			  lincelaux=int(auxl)
			end if
		  end if
		  !   (coluna)
		  if (colaux<pcol) then
			colcelaux=1
		  else
			auxc=colaux/real(pcol)
			auxc2=auxc-int(auxc)
			if (auxc2>0) then
			  colcelaux=int(auxc)+1
			else
			  colcelaux=int(auxc)
			end if
		  end if
		  if ((lincelaux==lincel).AND.(colcelaux==colcel)) then
			montante=0
			contacam=contacam+1     !pixel dentro da celula -> continuar caminho
			linaux2=linaux
			colaux2=colaux
			MaiorArea2=MaiorArea
		  else
			montante=1              !pixel fora da celula -> caminho encerrado
		  end if
		  if (contacam>LimiteCam) then
			montante=1              !jah atendeu o criterio caminho -> pode terminar caminho
		  end if
		end do     !fim do while montante=0

		if (contacam>LimiteCam) then
		  linpixex(lincel,colcel)=linexu
		  colpixex(lincel,colcel)=colexu
		  passo2=1      !pixel escolhido pela area atendeu ao criterio do caminho
		else
		  !pixel escolhido pela area NAO atendeu ao criterio do caminho
          
		  !vai ser analisado o HISTOGRAMA -> ver se o pixel testado eh o que drena o maior
		  !numero de pixels da celula (sim -> eh aceito; nao -> eh rejeitado)

		  	!loop para seguir caminho do fluxo de todos os pixels da celula
			do colpix=colpixini,colpixfin
			  do linpix=linpixini,linpixfin
				linaux=linpix
				colaux=colpix
				caminho=0
				!seguir caminho para jusante ateh sair da celula
				do while (caminho==0)
				  dd=dirpix(linaux,colaux)
				  linaux=linaux+dlin(dd)
				  colaux=colaux+dcol(dd)

				  !verifica se saiu area de estudo
				  if ((linaux<1).OR.(linaux>nlinpix).OR.(colaux<1).OR.(colaux>ncolpix)) then
					caminho=1
				  else

					  !verifica a qual celula o pixel pertence
					  !    (linha)
					  if (linaux<=plin) then
  						lincelaux=1
					  else
						auxl=linaux/real(plin)
						auxl2=auxl-int(auxl)
						if (auxl2>0) then
						  lincelaux=int(auxl)+1
						else
						  lincelaux=int(auxl)
						end if
					  end if
					  !   (coluna)
					  if (colaux<pcol) then
						colcelaux=1
					  else
						auxc=colaux/real(pcol)
						auxc2=auxc-int(auxc)
						if (auxc2>0) then
						  colcelaux=int(auxc)+1
						else
						  colcelaux=int(auxc)
						end if
					  end if
					  !checa se caminho saiu da celula
					  if ((lincelaux==lincel).AND.(colcelaux==colcel)) then
						!pixel estah dentro da celula -> caminho continua
						caminho=0
						!se o pixel estiver na borda, conta a passagem por ele
						if ((linaux==linpixini).OR.(linaux==linpixfin).OR.(colaux==colpixini).OR.(colaux==colpixfin)) then
						  !pixel na borda
						  histo(linaux,colaux)=histo(linaux,colaux)+1
						end if
					  else
						caminho=1
					  end if

				  end if !(fim if fora da area)

				end do !(fim do while caminho)
			  end do
			end do
    
			!loop pelos pixels da celula para ver qual drena a maior parte
			histomax=0
			do colpix=colpixini,colpixfin
			  do linpix=linpixini,linpixfin
				if (histo(linpix,colpix)>histomax) then
				  histomax=histo(linpix,colpix)
				  linexuH=linpix
				  colexuH=colpix
				end if
			  end do
			end do
    
	        !checa se o pixel testado para exutorio eh o escolhido pelo histograma
			!sim -> vai ser o pixel exutorio
			!nao -> vai ser rejeitado
			if ((linexu==linexuH).AND.(colexu==colexuH)) then
			  linpixex(lincel,colcel)=linexu
			  colpixex(lincel,colcel)=colexu
			  histocel(lincel,colcel)=histomax
			  passo2=1
			else
			  passo2=0
 				do colpix=colpixini,colpixfin
				  do linpix=linpixini,linpixfin
					histo(linpix,colpix)=0
				  end do
				end do
			end if



		end if !(fim if contacam)

  end do  !fim do while passo2=0

  end if !(fim if marcel=0)

     
  ! move a janela formada por plin x pcol pixels
  ! (vai para a proxima celula)	 
  if (colpixfin==ncolpix) then
    linpixini=linpixfin+1
	linpixfin=linpixfin+plin
	lincel=lincel+1
	colcel=1
    colpixini=1
	colpixfin=pcol
  else
    colcel=colcel+1
    colpixini=colpixfin+1
    colpixfin=colpixfin+pcol
  end if
  if (linpixfin>nlinpix) then
    fim=1
  else
    fim=0
  end if
end do  ! fim do while fim=0

!fim da identificao dos pixels exutorio
!!!write(*,*) "ok"
!!!write(*,*)


!>>>>>>>>>>>>>>>>>>>>>>>>> DETERMINACAO DAS DIRECOES DE FLUXO >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!!write(*,*) "3. definindo direcao de fluxo em cada celula..."

!celulas na borda drenam para fora da area de estudo e celulas vizinhas ao mar drenam para ele
!!!write(*,*) "<A.i>-atribuindo direcao pixels na borda ou vizinhos ao mar"
do colcel=1,ncolcel
  dircel(1,colcel)=H
  dircel(nlincel,colcel)=D
  do lincel=1,nlincel
    dircel(lincel,1)=F
	dircel(lincel,ncolcel)=B
    do i=1,8
 	  dd=ddaux(i)
	  laux=lincel+dlin(dd)
	  caux=colcel+dcol(dd)
	  if ((laux>=1).AND.(laux<=nlincel).AND.(caux>=1).AND.(caux<=ncolcel)) then
        if (marcel(laux,caux)==1) then  !o vizinho eh mar
		  vizmarcel(lincel,colcel)=1
	      dircel(lincel,colcel)=dd
	    end if !(fim if mar)
	  end if
	end do
  end do
end do
!!!write(*,*) "ok"
!!!write(*,*)

!loop para seguir caminho do fluxo
!
do lincel=1,nlincel
  do colcel=1,ncolcel
    write(*,*) lincel,colcel
    if ((marcel(lincel,colcel)==0).AND.(vizmarcel(lincel,colcel)==0)) then
		caminho=0
		linexu=linpixex(lincel,colcel)
		colexu=colpixex(lincel,colcel)
		dd=dirpix(linexu,colexu)
		!checa se a celula drena diretamente para fora da regiao estudada e atribui direcao 
		! para fora da regiao
		!(checa a posicao e direcao do pixel exutorio)
		if (linexu==nlinpix) then 
		  if (colexu==ncolpix) then
			dircel(lincel,colcel)=B
			caminho=1
		  else
			if (colexu==1) then
			  dircel(lincel,colcel)=F
			  caminho=1	
			else
			  if ((dd==C).OR.(dd==D).OR.(dd==E).OR.(dd==H)) then
				dircel(lincel,colcel)=D
	  			caminho=1
			  end if
			end if
		  end if
		end if
		if (linexu==1) then 
		  if (colexu==ncolpix) then
			dircel(lincel,colcel)=B
		  else  	  
			if (colexu==1) then
			  dircel(lincel,colcel)=F
			  caminho=1				
			else
			  if ((dd==A).OR.(dd==D).OR.(dd==G).OR.(dd==H)) then
				dircel(lincel,colcel)=H
				caminho=1
			  end if
			end if
		  end if
		end if
		if (colexu==ncolpix) then 
		  if ((linexu/=nlinpix).AND.(linexu/=1)) then	    
			if ((dd==A).OR.(dd==B).OR.(dd==C).OR.(dd==F)) then
			  dircel(lincel,colcel)=B
			  caminho=1
			end if
		  end if
		end if
		if (colexu==1) then 
		  if ((linexu/=nlinpix).AND.(linexu/=1)) then
			if ((dd==B).OR.(dd==E).OR.(dd==F).OR.(dd==G)) then
			  dircel(lincel,colcel)=F
			  caminho=1
			end if
		  end if
		end if


		exut=0
		linaux=linexu
		colaux=colexu

		! celula nao drena diretamente para fora da regiao estudada e
		! a busca vai ser iniciada ateh que [(area drenada pelo pixel exutorio encontrado) - (area 
		! drenada pixel exutorio cel inicial)] seja maior que um valor especificado
		do while (caminho==0)
		  diraux=dirpix(linaux,colaux)
		  foraviz=0
		  !determina o proximo pixel a jusante
		  do i=1,8
			dd=ddaux(i)
			if (diraux==dd) then
			  linaux=linaux+dlin(dd)
			  colaux=colaux+dcol(dd)
			end if
		  end do	
  
		 !identifica a qual celula o pixel encontrado pertence
		 !    (linha)
		 if (linaux<=plin) then
		   lincelaux=1
		 else
		   auxl=linaux/real(plin)
		   auxl2=auxl-int(auxl)
		   if (auxl2>0) then
			 lincelaux=int(auxl)+1
		   else
			 lincelaux=int(auxl)
		   end if
		 end if
		 !   (coluna)
		 if (colaux<pcol) then
		   colcelaux=1
		 else
		   auxc=colaux/real(pcol)
		   auxc2=auxc-int(auxc)
		   if (auxc2>0) then
			 colcelaux=int(auxc)+1
		   else
			 colcelaux=int(auxc)
		   end if
		 end if
		 
    
		 !verifica se saiu para alem das 8 celulas vizinhas
		 ! (entao a direcao nao serah alterada, permanecendo (i) a referente ao ultimo 
		 ! pixel exutorio encontrado, caso tenha sido encontrado algum ou (ii) a direcao
		 ! relativa a celula pela qual o caminho saiu da vizinhanca) 
		 if ((lincelaux>lincel+1).OR.(colcelaux>colcel+1)) then
		   foraviz=1
		   caminho=1
		 end if
		 if ((lincelaux<lincel-1).OR.(colcelaux<colcel-1)) then
		   foraviz=1
		   caminho=1
		 end if
		 ! caso em que saiu sem passar por nenhum pixel exutorio
		 if ((foraviz==1).AND.(exut==0)) then  
		       
			!direcao atribuida em funcao de por onde o caminho saiu da viz
			 if (lincelaux>lincel+1) then  !(saiu por baixo da viz...)
			   if (colcelaux<=colcel-1) then !(e esquerda)
				 dircel(lincel,colcel)=E
			   end if
			   if (colcelaux==colcel) then   !(e centro)
				 dircel(lincel,colcel)=D
			   end if
			   if (colcelaux>=colcel+1) then !(e direita)
				 dircel(lincel,colcel)=C
			   end if
			 end if
			 if (lincelaux<lincel-1) then !(saiu por cima da viz..)
			   if (colcelaux<=colcel-1) then !(e esquerda)
				 dircel(lincel,colcel)=G
			   end if
			   if (colcelaux==colcel) then !(e centro)
				 dircel(lincel,colcel)=H
			   end if
			   if (colcelaux>=colcel+1) then !(e direita)
				 dircel(lincel,colcel)=A
			   end if
			 end if
			 if (colcelaux>colcel+1) then !(saiu pela direita da viz..)
			   if (lincelaux<=lincel-1) then !(e acima)
				 dircel(lincel,colcel)=A
			   end if
			   if (lincelaux==lincel) then  !(e centro)
				 dircel(lincel,colcel)=B
			   end if
			   if (lincelaux>=lincel+1) then !(e abaixo)
				 dircel(lincel,colcel)=C
			   end if
			 end if
			 if (colcelaux<colcel-1) then !(saiu pela esquerda da viz..)
			   if (lincelaux<=lincel-1) then
				 dircel(lincel,colcel)=G  !(e acima)
			   end if
			   if (lincelaux==lincel) then
				 dircel(lincel,colcel)=F !(e centro)
			   end if
			   if (lincelaux>=lincel+1) then
				 dircel(lincel,colcel)=E  !(e abaixo)
			   end if
			 end if
		  end if
		  ! caso em que saiu MAS passou por algum pixel exutorio
		  if ((foraviz==1).AND.(exut==1)) then    
			dircel(lincel,colcel)=dircel(lincel,colcel)
		  end if
     

		 !verifica se saiu para fora da regiao de estudo
		 if ((linaux<1).OR.(linaux>nlinpix)) then
		   foraviz=1
		   caminho=1
		   if (exut==0) then   !significa que nao passou por nenhum pix exutorio
			 dircel(lincel,colcel)=dirpix(linexu,colexu)
		   else
			 dircel(lincel,colcel)=dircel(lincel,colcel)
		   end if
		 end if 
		 if ((colaux<1).OR.(colaux>ncolpix)) then
		   foraviz=1
		   caminho=1
		   if (exut==0) then
			 dircel(lincel,colcel)=dirpix(linexu,colexu)
		   else
			 dircel(lincel,colcel)=dircel(lincel,colcel)
		   end if
		 end if 



		 !caso nao tenha saido da vizinhanca nem da regiao de estudo, verifica se o pixel eh 
		 !exutorio e checa criterio da area, atribuindo a direcao em caso de atendimento
		 if (foraviz==0)then
		   !verifica se o pixel eh exutorio da celula
  		   if (linaux==linpixex(lincelaux,colcelaux)) then
			 if (colaux==colpixex(lincelaux,colcelaux)) then  !eh exutorio!
				 exut=1
			   !verifica area de drenagem incremental
				 AreaInc=aacpix(linaux,colaux)-aacpix(linexu,colexu)		    
				 if (AreaInc>LimiteInc) then
				   caminho=1
				 end if
				 deltacol=colcelaux-colcel
				 deltalin=lincelaux-lincel
				 if (deltacol==1) then
				   if (deltalin==-1) then
					 dircel(lincel,colcel)=A
				   end if
				   if (deltalin==0) then
					 dircel(lincel,colcel)=B
				   end if
				   if (deltalin==1) then
					 dircel(lincel,colcel)=C
				   end if
				 end if
				 if (deltacol==0) then
				   if (deltalin==-1) then
					 dircel(lincel,colcel)=H
				   end if
				   if (deltalin==0) then
					 dircel(lincel,colcel)=0
				   end if
				   if (deltalin==1) then
					 dircel(lincel,colcel)=D
				   end if
				 end if
				 if (deltacol==-1) then
				   if (deltalin==-1) then
					 dircel(lincel,colcel)=G
				   end if
				   if (deltalin==0) then
					 dircel(lincel,colcel)=F
				   end if
				   if (deltalin==1) then
					 dircel(lincel,colcel)=E
				   end if
				 end if             
			 end if
		   end if
		 end if
		end do
	end if !(fim if marcel and vizmarcel)
  end do
end do

!!!write(*,*) "ok"
!!!write(*,*)



!>>>>>>>>>>>>>>>>>>> VERIFICACAO E CORRECAO DE CRUZAMENTOS NAS DIRECOES DE FLUXO >>>>>>>>>>>>

write(*,*) "4. verificacao e correcao de 'cruzamentos'..."

!janela 2 x 2 celulas
!    | dircelSE  dircelSD |
!    | dircelIE  dircelID |
!     
cruzamento=0
correcao=0
do while (correcao==0)   !soh termina quando nao houver nenhuma correcao no mesmo loop
  do lincel=1,nlincel-1   
    do colcel=1,ncolcel-1
      if ((marcel(lincel,colcel)==0).AND.(vizmarcel(lincel,colcel)==0)) then

		  dircelSE=dircel(lincel,colcel)
		  dircelSD=dircel(lincel,colcel+1)
		  dircelIE=dircel(lincel+1,colcel)
		  dircelID=dircel(lincel+1,colcel+1)
		  !armazena areas drenagem dos pixels exutorios de cada celula da janela
		  aacpixSE=aacpix(linpixex(lincel,colcel),colpixex(lincel,colcel))
		  aacpixSD=aacpix(linpixex(lincel,colcel+1),colpixex(lincel,colcel+1))
		  aacpixIE=aacpix(linpixex(lincel+1,colcel),colpixex(lincel+1,colcel))
		  aacpixID=aacpix(linpixex(lincel+1,colcel+1),colpixex(lincel+1,colcel+1))
		  if ((dircelIE==A).AND.(dircelID==G)) then
			cruzamento=cruzamento+1
			if ((dircelSE/=D).AND.(aacpixIE<aacpixID)) then
			  dircelIE=H
			  dircel(lincel+1,colcel)=dircelIE
			else
			  dircelID=H
			  dircel(lincel+1,colcel+1)=dircelID
			end if
		  end if
		  if ((dircelSE==C).AND.(dircelSD==E)) then
			cruzamento=cruzamento+1
			if ((dircelIE/=H).AND.(aacpixSE<aacpixSD)) then
			  dircelSE=D
			  dircel(lincel,colcel)=dircelSE
			else
			  dircelSD=D
			  dircel(lincel,colcel+1)=dircelSD
			end if
		  end if
		  if ((dircelSD==E).AND.(dircelID==G)) then
			cruzamento=cruzamento+1
			if ((dircelSE/=B).AND.(aacpixSD<aacpixID)) then
			  dircelSD=F
			  dircel(lincel,colcel+1)=dircelSD
			else
			  dircelID=F
			  dircel(lincel+1,colcel+1)=dircelID
			end if
		  end if
		  if ((dircelSE==C).AND.(dircelIE==A)) then
			cruzamento=cruzamento+1
			if ((dircelSD/=F).AND.(aacpixSE<aacpixIE)) then
			  dircelSE=B
			  dircel(lincel,colcel)=dircelSE
			else
			  dircelIE=B
			  dircel(lincel+1,colcel)=dircelIE
			end if
		  end if
		 end if !(fim if marcel and vizmarcel)
		end do
	   end do
	  write(*,*) cruzamento
	  if (cruzamento==0) then
		correcao=1
	  else
		correcao=0
	  end if
	  cruzamento=0
end do  !fim do while correcao=0

!!!write(*,*) "ok"
!!!write(*,*)


!>>>>>>>>>>>>>>>>>>>> analise das direcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!!!write(*,*) "4.2. analisando situacoes especificas..."

cont=0
marcacel=0
do colcel=1,ncolcel
  do lincel=1,nlincel
    if (marcel(lincel,colcel)==0) then
	  linexu=linpixex(lincel,colcel)
	  colexu=colpixex(lincel,colcel)
	  AreaMont=aacpix(linexu,colexu)
	  dd=dircel(lincel,colcel)
	  lincelaux=lincel+dlin(dd)
	  colcelaux=colcel+dcol(dd)
	  if ((lincelaux>=1).AND.(lincelaux<=nlincel).AND.(colcelaux>=1).AND.(colcelaux<=ncolcel)) then
	  	if (marcel(lincelaux,colcelaux)==0) then
			linaux=linpixex(lincelaux,colcelaux)
			colaux=colpixex(lincelaux,colcelaux)
			AreaJus=aacpix(linaux,colaux)
			if (AreaJus<AreaMont) then
			  cont=cont+1
			  marcacel(lincel,colcel)=1
			end if
		end if
	  end if
	end if
  end do
end do

!!!write(*,*) "numero de situacoes onde Ajus<Amont: ",cont
!!!write(*,*)



!>>>>>>>>>>>>>>>>>>>> PROCESSAMENTO DE ARQUIVOS DE SAIDA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!!!write(*,*) "5. escrevendo arquivo de saida..."
!!!write(*,*)
!!!write(*,*) "DirBAIXA.rst - (int/bin) - direcoes de fluxo de baixa resolucao"
!!!write(*,*) "PixelExu.dat - (ascii) - matriz com linha e coluna dos pixels exutorios"
!!!write(*,*) "PixelExu.RST - (int/bin) - imagem ALTA resolucao indicando pixels exutorios"


    
! escrita de arquivo de direcoes - arquivo binario/inteiro
!!!open(60,file='DIRBAIXA.rst',status='unknown',form='unformatted',access='direct',RECL=2*ncolcel)
!!!do lincel=1,nlincel
!!!  write(60,REC=lincel) (dircel(lincel,colcel),colcel=1,ncolcel)
!!!end do
!!!close(60)

!!!open(42,file='pixelexu.dat')
!!!do lincel=1,nlincel
!!!  write(42,'(2000i6)') (linpixex(lincel,colcel),colpixex(lincel,colcel),colcel=1,ncolcel)
!!!end do

! geracao de matriz onde pixels exutorios tem 1 e demais tem 0
exu=0
do lincel=1,nlincel
  do colcel=1,ncolcel
    linexu=linpixex(lincel,colcel)
	colexu=colpixex(lincel,colcel)
	exu(linexu,colexu)=1
  end do
end do

!escrita de arquivo para visualizar pixels exutorios
!!!open(69,file='pixelexu.rst',status='unknown',form='unformatted',access='direct',RECL=2*ncolpix)
!!!do linpix=1,nlinpix
!!!  write(69,REC=linpix) (exu(linpix,colpix),colpix=1,ncolpix)
!!!end do
!!!close(69)


!!!write(*,*) "ok"
!!!write(*,*)

!!!write(*,*) "Programa DirFluxo terminou! (tecle enter)"
!!!write(*,*)



! Copy OUTPUT data - for R compatibility
do j = 1, nlinpix
   do i = 1, ncolpix
      outlet(i,j) = exu(j,i)
   end do
end do
do j = 1, nlincel
   do i = 1, ncolcel
      dirl(i,j) = dircel(j,i)
   end do
end do
do j = 1, nlincel
   do i = 1, ncolcel
      outrow(i,j) = linpixex(j,i)
   end do
end do
do j = 1, nlincel
   do i = 1, ncolcel
      outcol(i,j) = colpixex(j,i)
   end do
end do

end subroutine fluxdir

!!!end Program DirFluxo4

