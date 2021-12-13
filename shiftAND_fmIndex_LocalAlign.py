from collections import namedtuple
from bitarray import bitarray
import time
import textwrap


class FMNucleo:
    def __init__(self, base, index):
        self.base = base
        self.index = index


class FatsaFile:
    # başlık satırının uzunluğu, indexleme yaparken bulunan index + header_linelength
    # ile doğru index numarasına erişiriz
    header_linelength = None

    def __init__(self, filename):
        self.filename = filename

    # fatsa dosyası içerisindeki metnin tek satır hali
    # ve dosyanın ilk satırının uzunluğunu hesaplarız
    def readfile(self):
        file = open(self.filename, 'r')
        self.header_linelength = len(next(file))
        text = ""
        for line in file:
            for char in line.strip():
                text += char
        return text


class Operations:
    # ön işlemler ve arama işlemlerinin sürelerini tutmak için
    # oluşturulan değişkenler
    shiftand_preprocesstime = 0
    shiftand_totalsearchtime = 0
    fm_preprocesstime = 0
    fm_totalsearchtime = 0

    # smith waterman patterını
    sw_pattern = ""

    def __init__(self, patterns, text, header_linelength):
        self.patterns = patterns
        self.text = text
        self.header_linelength = header_linelength

    # shift and algoritması
    # öncelikli olarak verilen tüm patternler için alfabe bilgisi elde edilir
    # sonrası bu alfabe bilgisine göre metin içerisinde aramalara başlanır
    # okunan her bir karakterin örüntü içerisindeki bit dizisi karşılığı elde edilir
    # shifter değişkeni her işlemin sonundaki durumu ifade eden bit dizisidir
    # shifter[-1] == 0 olduğunda örüntü metin içerisinde yakalanmış demektir
    def shiftand(self):
        # ön işlem süresi hesaplanır
        # alphabet_and_nucleobits her bir pattern alfabesine karşılık gelen karakterlerin
        # pattern içerisindeki bit dizisini ttuar
        preprocess_start = time.time()
        alphabet_and_nucleobits = self.shiftand_getalphabet_and_nucleobits()
        preprocess_end = time.time()
        self.shiftand_preprocesstime = preprocess_end - preprocess_start

        # bitarray olarak oluşturulan shifterın tüm bitleri
        # önce sıfıra eşitlenir. sonra ilk bit 1 yapılır
        # matches ile örüntünün kaç kez bulunduğu bilgisini tutarız
        # shifter_1 değişkeni shifter'ın her işlem sonunda 1 sağa
        # kaydırılmasının ardından ilk bitinin 1 olması için üretilmiştir
        # 10000.. şeklinde başlayan bu bit dizisi ile OR yapılan shifter
        # sonucunda 1 ile başlamış olur
        # patternIndex ile okunan örüntünün index numarası alınır
        # bu numara ile alphabet_..._bits içerisinde okunan karakterin
        # bit dizisi aranır. eğer bulunursa shifter ile AND işlemine tabi
        # tutulur. tüm karakterler için bu AND işlemi yapılırken eğer
        # shifter[-1] == 0 örüntü bulunmuş olur
        for pattern in self.patterns:
            preprocess_start = time.time()
            shifter = bitarray(len(pattern))
            shifter.setall(0)
            shifter[0] = 1
            shifter_1 = shifter.copy()
            patternIndex = self.patterns.index(pattern)
            matches = 0
            preprocess_end = time.time()
            self.shiftand_preprocesstime += preprocess_start - preprocess_end

            for textindex, t in enumerate(self.text, 0):
                pattern_searchstart = time.time()
                nucleo = next((nucleo for nucleo in alphabet_and_nucleobits[patternIndex] if nucleo.base == t), None)
                if nucleo is not None:
                    shifter = shifter.__and__(nucleo.bits)
                    if shifter[-1] == 1:
                        pattern_searchend = time.time()
                        print(pattern, self.header_linelength + textindex - len(pattern), ". ve", self.header_linelength + textindex, ". karakterler arasında bulundu")
                        self.shiftand_totalsearchtime += pattern_searchend - pattern_searchstart  # her bir başarılı aramanın bulma süresi
                        matches += 1
                        shifter.setall(0)
                        shifter[0] = 1
                    shifter = shifter >> 1
                    # python sağa kaydırmada birinci sıradan 1 göndermediği için
                    # ilk elemanı 1 olan (geri kalanlar 0) sabit shifter_1 ile shifter orlandığında
                    # shifterın başına 1 eklenmiş olur
                    shifter = shifter.__or__(shifter_1)
                else:
                    # eğer okunan karakter alfabe içerisinde yoksa
                    # o anki text index sırasından sonrası için shifter sıfırlanır
                    # çünkü patternın bulunması için okunan o anki karakterin
                    # pattern içerisinde olması gerekmektedir
                    shifter.setall(0)
                    shifter[0] = 1
            print("-------------------------------------------------------------------------")
            print(pattern, "dosya içerisinde toplam ", matches, " kez bulundu")
            print(pattern, "için yaklaşık toplam arama süresi =", self.shiftand_totalsearchtime)
        print("-------------------------------------------------------------------------")

    # shift and işlemi için verilen örüntülerin alfabesini ve
    # bu alfabelerdeki karakterlerin örüntüye içinde tekrar ettiği
    # yerlere karşılık gelen 1 ve etmeyen 0'ları üreten
    # ve bu tüm alfabe bit dizilerini nükleobaazlarıyla birlikte
    # tuple listesi olarak döndüren fonksiyon
    def shiftand_getalphabet_and_nucleobits(self):
        alphabet_nucleos = []
        Pair = namedtuple("nucleo", ["base", "bits"])
        for patterns in self.patterns:
            # her seferinde alfabenin aynı sırada oluşması için sort işlemi
            alphabet = sorted(set(patterns))
            nucleos = []
            for base in alphabet:
                nucleobase_bits = bitarray()
                for patternbase in patterns:
                    if base == patternbase:
                        nucleobase_bits.append(1)
                    else:
                        nucleobase_bits.append(0)
                nucleos.append(Pair(base, nucleobase_bits))
            alphabet_nucleos.append(nucleos)
        return alphabet_nucleos

    # fm_index fonksiyonu
    # fm index için ilk olarak metnin rotasyon tablosu, burrows wheeler transormasyonu elde edilir
    # bu transformasyonun ilk satırı kullanılarak C = c_n_occurrences değerleri elde edilir
    # C burada örneğin $AAACCCCGGGTTTT şeklinde sıralı olan ilk BWT satırı için
    # her bir alfabe karakterinin ilk hangi indexte görüldüğünün bilgisni tutar
    # üst örnek için C = [0, 1, 4, 8, 11] olarak (alfabe sırasına göre) elde edilir
    # sonrasında eğer örüntü bulunduğunda metin içerisinde hangi indexte olduğuyla ilgileniyorsak
    # bir suffix_array üretilir. suffix array BWT nin son satırına göre alfabetik olarak sıralanıp
    # her bir satır içerisinde $ işaretinden sonraki kısımlar, yani son eklerin bilgisini tutar
    # bu sayede yapılan sonraki işlemlerden elde edilen aralık değerleri ile bulunan örüntülerin
    # hangi indexlerde olduğunu hesaplayabiliriz. suffix arrayde aynı şekilde başlayan alt metinler
    # yan yana bulunacağı için elde edilen aralık değeri kullanılarak suffix_array'de karşılık gelen
    # indexlerdeki sayıların her biri örüntünün metin içerisindeki index bilgisini verir
    # yine BWT'nin son sütunu kullanılarak bir Occ matrisi elde edilir
    # Occ matrisi son satırın içerisinde yukarıdan aşağıda her bir karakterin kaç kez görüldüğü
    # bilgisini tutar. Bu bilgi her bir index için kayıt edilir.
    # 2 karakter ve $ işaretinden oluşan bir alfabenin Occ si [0,0,1], [0,1,1], [0,2,1] [0, 3, 1] ... gibi bir
    # şekilde olacaktır. bu bilgi bizlere o anki index için karakterlerin o ana kadar kaç kez okunduğu
    # bilgisini verir
    # algoritma geliştirilirken çeşitli araştırma makalelerinden faydalanılmıştır
    # OffScan: a universal and fast CRISPR offtarget sites detection tool
    # Yingbo Cui1, Xiangke Liao1, Shaoliang Peng2,3, Tao Tang1, Chun Huang1* and Canqun Yang1
    # çalışma içerisinde fm index algoritmasına yer verildiği için implemantasyon bu örnek
    # üzerinden denenmiştir

    # sp değişkeni sondan okuması yapılan örüntünün o ana kadar okunan karakterlerinin
    # metin içerisindeki aralığının başlangıç konumunu ifade eder.
    # ep ise bitiş konumunu ifade eder
    # sonuçta yaptığımız sıralamalardan dolayı birbirine benzer olan yapılar yan yana kümelendiği için
    # başlangıç noktası bitişten büyük olana kadar (bu durumda yalnızca 1 tane eşleşme olurdu ancak başka
    # bir değer daha kullanıyoruz) ve sondan okuması yapılan örüntünün ilk karakterine ulaşıncaya kadar
    # bu işlem devam eder. ne zaman ki son karakter okunur ve sp <= ep dir, o anda artık
    # örüntünün metin içerisindeki tüm tekrarlarının toplam sayısını elde edebiliriz
    # ep - sp + 1
    # ep ve sp değerlerini kullanarak suffix_array yardımıyla örüntülerin metin içerisindeki
    # index numaralarını da bulabiliriz. yalnızca tekrar sayısı ile ilgileniyorsak vr buna ihtiyacımız
    # yoksa ise ön işlem süresinde bir kazanç elde ederiz.
    def fm_index(self):
        preprocess_start = time.time()

        self.text += "$"
        transformlist = []
        for i in range(0, len(self.text)):
            transformlist.append(self.text[i:] + self.text[0:i])
        sortedrotations = sorted(transformlist)

        alphabet = ["$", "A", "C", "G", "T"]

        # eğer text içerisinde patternlerin bulunduğu konumu belirlemek istiyorsak suffix_array oluşturmalıyız
        # bu işlemin yapılması zamanı fm_index çerçevesinde fazlaca arttıracaktır
        suffix_array = self.fm_getsuffixarray(sortedrotations)
        bwt_lastcol = self.fm_getbwt(sortedrotations)
        c_noccurrences = self.fm_get_cnoccurrences(sortedrotations, alphabet)
        occ = self.fm_get_occ(bwt_lastcol)

        preprocess_end = time.time()
        self.fm_preprocesstime = preprocess_end - preprocess_start

        # sp suffix array içerisinde bulunan patternın başlangıç noktası
        # ep suffix array içerisinde bulunan patternın bitiş noktası
        # ep ve sp Occ matrisi ve C listesi kullanıarak hesaplanır
        # Occ hangi indexte kaç kez karakter tekrar olduğu, C karakterlerin sıralı yapıdaki (bwt ilk kolonun sıralı
        # hali) ilk okunma index,
        for pattern in self.patterns:
            fm_searcstart = time.time()
            i = len(pattern) - 1
            c_index = alphabet.index(pattern[i])
            sp = c_noccurrences[c_index]
            if c_index == len(alphabet) - 1:
                ep = sp + occ[-1][-1] - 1
            else:
                ep = c_noccurrences[c_index + 1]
            while sp <= ep and i >= 0:
                if sp <= ep and i == 0:
                    fm_searchend = time.time()
                    print(pattern, "dosya içerisinde",  ep - sp + 1, "kez bulundu")
                    print(pattern, "için yaklaşık toplam arama süresi =", fm_searchend - fm_searcstart)
                    print("--------------------------------------------------------------")
                    self.fm_totalsearchtime += fm_searchend - fm_searcstart
                    # aranan patterin text içerisinde hangi konumda olduğunu görmek istiyorsak
                    # if ep - sp != 0:
                    #     for sa_index in range(sp, ep, 1):
                    #         print(pattern, suffix_array[sa_index] + self.header_linelength - 1, ". ve ", suffix_array[sa_index] + len(pattern) + self.header_linelength - 1, ". karakterler arasında bulundu")
                    # else:
                    #     for sa_index in range(sp, ep+1, 1):
                    #         print(pattern, suffix_array[sa_index] + self.header_linelength - 1, ". ve ", suffix_array[sa_index] + len(pattern) + self.header_linelength - 1, ". karakterler arasında bulundu")
                kk = alphabet.index(pattern[i - 1])
                sp = c_noccurrences[kk] + occ[sp - 1][kk]
                ep = c_noccurrences[kk] + occ[ep][kk] - 1
                i -= 1

    # herhangi verilen bwt için suffix array oluşturur
    # suffix_array listesine her bir $ işaretinden sonraki
    # kısımları ekler
    @staticmethod
    def fm_getsuffixarray(sortedrotations):
        suffix_array = []
        for rotation in sortedrotations:
            suffix_array.append(len(rotation) - rotation.find("$") - 1)
        return suffix_array

    # bwt'nin son sütununu döndürür
    @staticmethod
    def fm_getbwt(sortedrotations):
        bwt_lastcol = ""
        for rotation in sortedrotations:
            bwt_lastcol += rotation[-1]
        return bwt_lastcol

    # bwt sıralanmış ilk sütunu içerisindeki karakterlerin
    # ilk hangi indexte görüldüğü bilgisini tutan
    # C (c_n_occurrences) bilgisini döndürür
    @staticmethod
    def fm_get_cnoccurrences(sortedrotations, alphabet):
        bwt_firstcol = ""

        for rotation in sortedrotations:
            bwt_firstcol += rotation[0]

        c_noccurrences = []
        for letter in alphabet:
            c_noccurrences.append(bwt_firstcol.find(letter))

        return c_noccurrences

    # bwt son sütununu kullanarak sütun içerisinde
    # yukarından aşağı okuma yapılırken görülen her bir
    # karakterin o ana kadar kaç kez görüldüğünü ve önceki durumları da
    # liste şeklinde tutan Occ matrisini döndürür
    @staticmethod
    def fm_get_occ(bwt_lastcol):
        c_a = 0  # indexe göre A görülme sayısı
        c_c = 0  # indexe göre C görülme sayısı
        c_t = 0  # indexe göre T görülme sayısı
        c_g = 0  # indexe göre G görülme sayısı
        c_d = 0  # indexe göre $ görülme sayısı
        occ = []
        for i in range(0, len(bwt_lastcol)):
            if i == 0:
                if bwt_lastcol[i] == "$":
                    c_d += 1
                elif bwt_lastcol[i] == "A":
                    c_a += 1
                elif bwt_lastcol[i] == "C":
                    c_c += 1
                elif bwt_lastcol[i] == "G":
                    c_g += 1
                else:
                    c_t += 1
                occ.append([c_d, c_a, c_c, c_g, c_t])
            else:
                if bwt_lastcol[i] == "$":
                    c_d += 1
                elif bwt_lastcol[i] == "A":
                    c_a += 1
                elif bwt_lastcol[i] == "C":
                    c_c += 1
                elif bwt_lastcol[i] == "G":
                    c_g += 1
                else:
                    c_t += 1
                occ.append([c_d, c_a, c_c, c_g, c_t])
        return occ

    # smith waterman local alignment işlemlerini yapan fonksioyon
    # öncelikle (len(metin) + 1 * len(örüntü) + 1) i,j elemana sahip bir
    # sıfırlardan oluşan rank matrisi üretilir. global allignmenttan farklı
    # olarak bu matriste negatif değerler olmaz
    # matris içerisinde; benzerlik yolu bilgisini tutar
    # tersten okunması yapılan örüntü ile denk gelen karakterlerin
    # olduğu yerlerde diyagonal, olmadığı yerlerde duruma göre
    # yukarı veya sola gidiyorsa metine ya da örüntüye - işareti
    # koyarak benzerliği bulan fonksiyon
    # ilk başta üretilen matrisin en büyük değere (best_rank) sahip
    # hücresinin değeri başlangıç i,j konumunu ifade eder
    # best_i bu hücrenin satır bilgisi
    # best_j bu hücrenin sütun bilgisi
    # bu bilgilere geri okuma yaparken ihtiyacımız olacaktır.
    def sw_local(self, match, mismatch, gap):
        sub_texts = self.slice_text(4)
        scores = []
        for sub_text in sub_texts:
            for i in range(len(sub_text) + 1):
                scores.append([0] * (len(self.sw_pattern) + 1))

            best_rank = 0
            best_i = 0
            best_j = 0

            # rank matrisini oluşturan yapı
            # verilen match, mismatch ve gap değerlerine göre
            # geçişler diyagonalse ve match varsa, mismatch varsa
            # geçişler yatay yönlü ise + gap
            # geçişler dikey yönlü ise + gap
            # yaparak max(0, D, U, L) alınarak
            # bir sonraki hücrenin değeri hesaplanır
            # bu değer max içerisinde 0 olduğu için (lokal hizalama)
            # negatif olamaz
            for i in range(1, len(sub_text) + 1):
                for j in range(1, len(self.sw_pattern) + 1):
                    D = scores[i - 1][j - 1] + match if sub_text[i - 1] == self.sw_pattern[j - 1] else mismatch
                    U = scores[i - 1][j] + gap
                    L = scores[i][j - 1] + gap
                    scores[i][j] = max(0, D, U, L)
                    # en yüksek skora sahip hücre bulunur
                    # bu hücrenin index bilgisi tutulur
                    if scores[i][j] >= best_rank:
                        best_rank = scores[i][j]
                        best_i = i
                        best_j = j

            self.sw_backtrack(scores, sub_text, best_i, best_j, match, mismatch, gap)

    # rank matrisi üzerinde örüntünün son karakterinden okunarak bir geri yönlü
    # yol oluşturulması yapılır
    # bu sayede bükülme olan yerlere "-" gelecek şekilde
    # ve eşleşme olan yerler aynı kalacak şekilde
    # local alignment yapılmış olur
    def sw_backtrack(self, scores, text, reverse_i, reverse_j, match, mismatch, gap):
        local_text = ""
        local_pattern = ""
        while (reverse_i > 0 or reverse_j > 0) and scores[reverse_i][reverse_j] != 0:
            D = 0  # diyagonal geçiş
            U = 0  # dikey geçiş, up
            L = 0  # yatay geçiş, left

            # eğer score[0,0] a gelmediysek
            if reverse_i > 0 and reverse_j > 0:
                # okunan hücre diyagonal mi üretildi onu bulabilmek adına
                # eğer diyagonal olsaydı hangi değeri üretirdi onu D'de tutarız
                D = scores[reverse_i - 1][reverse_j - 1] + match if text[reverse_i - 1] == self.sw_pattern[reverse_j - 1] else mismatch
            if reverse_i > 0:
                # eğer i henüz bitmediyse dikey yönde gelmiş olabiliriz
                # U eğer dikey gelinmiş olsaydı gap olacağı için
                # hücrenin değeri bir üst hücre + gap ten mi oluştu
                # onu öğrenmek için U da bu değeri tutarız
                U = (scores[reverse_i - 1][reverse_j] + gap)
            if reverse_j > 0:
                # U da olduğu gibi bu hücre değerini soldan mı ürettik
                # onu anlayabilmek için bir sol hücre + gap kadar
                # bir değer belirleriz
                L = (scores[reverse_i][reverse_j - 1] + gap)
            # eğer hücre diyagonal oluştuysa doğal olarak D puanı U ve L den
            # daha yüksek olacaktır
            # ve eşleşme diyagonal oluştuğu için hem text hem de örüntü tutucu boş stringe
            # okunan karakteri ekleriz
            if D >= U and D >= L:
                local_text = text[reverse_i - 1] + local_text
                local_pattern = self.sw_pattern[reverse_j - 1] + local_pattern
                reverse_i -= 1
                reverse_j -= 1
            # eğer U en büyük puana sahipse hücreye bir üst hücreden gelmişiz demektir
            # ve ilk satır text bilgisi olduğu için o yönde bir bükülme olur
            # yani karakter textten okunur ancak örüntüden okunmaz onun yerine "-" konulur
            elif U >= L:
                local_text = text[reverse_i - 1] + local_text
                local_pattern = '-' + local_pattern
                reverse_i -= 1
            # son durum olarak L > ise karakteri soldan okumuşuz demektir ve sol kısım
            # örüntüyü ifade ettiği için textten denk gelen yere "-" koyarız
            else:
                local_pattern = self.sw_pattern[reverse_j - 1] + local_pattern
                local_text = '-' + local_text
                reverse_j -= 1
        # local_text text içerisinde benzerlik olan bölgedeki benzerliği ifade eder
        # local_pattern ise pattern içerisinde texte denk en yüksek benzerliğe sahip olan
        # yeri ifade eder
        print("text içi:", local_text)
        print("pattern :", local_pattern)
        print("---------------------------------------------------------------")

    # ödev tanımında stringin 4 parçaya bölünmesi istenmiştir
    # tüm textin uzunluğu 4 e bölünerek elde edilen parçalanmış textleri
    # döndürür. textwrap kütüphanesi ile gerçekleştirilmiştir
    # ilk parametre bölünecek olan stringin ismi, ikinci parametre bölünecek olan
    # parçanın uzunluğu ne kadar onu ifade eder
    # bizim durumumuzda len(self.text) + 2 / 4 - - - - şeklinde
    # parçalanmış textleri oluşturur
    # +2 5. bir parça oluşmaması için eklenmiştir
    # 2 den fazla herhangi bir i <= k sayısı bu durumu sağlar
    # k = texti 3 parçaya böldüren sayı
    def slice_text(self, split_into_n):
        sub_texts = textwrap.wrap(self.text, int((len(self.text) + 2) / 4))
        return sub_texts


