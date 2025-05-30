# Changelog

## [0.9.3](https://github.com/Danderson123/Amira/compare/v0.9.2...v0.9.3) (2025-05-29)


### Miscellaneous Chores

* version bump ([f4e38d8](https://github.com/Danderson123/Amira/commit/f4e38d897644b05680ac7e66088a23913bfb62f3))


### Documentation

* add citation ([736de7c](https://github.com/Danderson123/Amira/commit/736de7c690544a74a0bf4fc9e09f71249c0adbdc))
* add conda install commands ([c5a6e6e](https://github.com/Danderson123/Amira/commit/c5a6e6e82b00bc7cc5dd452eb3b059e0f0a664bf))
* update README ([69a9978](https://github.com/Danderson123/Amira/commit/69a99786fdf6fb7ccf9758d9a6af32cb7a6fe010))


### Bug Fixes

* add KP plasmid gene files ([4ca3672](https://github.com/Danderson123/Amira/commit/4ca3672a72acd73ad4d076ef6a2104dd6c012c14))
* add SA and SP plasmid gene files ([0b8683e](https://github.com/Danderson123/Amira/commit/0b8683ea35ffe48f73a4b2f9d9ddfe0520b6728c))
* fix for underestimation of relative read depth in some cases ([347e5d4](https://github.com/Danderson123/Amira/commit/347e5d444429b6427955e23ae96619c47b012135))
* improper filtering of low identity genes ([4a3ddf5](https://github.com/Danderson123/Amira/commit/4a3ddf5333e5d15d514a64690f0d8c90d89b3bd6))
* pre-commit formatting ([2c2950a](https://github.com/Danderson123/Amira/commit/2c2950a514662ce67172a000d2fe320b02bee700))
* properly report coverage across genes with no fully covering reads ([2585180](https://github.com/Danderson123/Amira/commit/25851806bd432ba195a0d500690b7e1075a1fa7a))
* remove debug line ([c67ec67](https://github.com/Danderson123/Amira/commit/c67ec67da1e3bac89ab41767ce6dc7fc79e9aa42))
* undo coverage check change ([74e85ba](https://github.com/Danderson123/Amira/commit/74e85bab37f49d0d8e651f10c685ea3e663093fd))
* update setuptools version ([716b405](https://github.com/Danderson123/Amira/commit/716b4051ba8a433d5e50cc87d9a125364e566f09))

## [0.9.2](https://github.com/Danderson123/Amira/compare/v0.9.1...v0.9.2) (2025-05-16)


### Bug Fixes

* correct shortest paths first ([eeb77b5](https://github.com/Danderson123/Amira/commit/eeb77b51c9a64f7cb5599c98cd1734f0547be1f4))
* missing arguements when genotyping promoters ([f47e5e8](https://github.com/Danderson123/Amira/commit/f47e5e8df1baaa92105a7178da18cf1402efb952))
* reduce sourmash k to marginally improve precision ([94fc4ef](https://github.com/Danderson123/Amira/commit/94fc4ef55b77c8d9df111c0e2a0d8192b3709618))
* remove overwriting of JSON annotations ([6e71c13](https://github.com/Danderson123/Amira/commit/6e71c1329e420fc937fb74509e5aeecc7a48b0a9))
* update unittests for new sourmash params ([b0d951e](https://github.com/Danderson123/Amira/commit/b0d951e10fe2a7f2f21d2d282b24cbdf7e44db0c))
* use a single core to build the graph in the graph correction stage ([06cf34d](https://github.com/Danderson123/Amira/commit/06cf34dce6f492dd74f7d35fbf36c5c49ec5cdf5))

## [0.9.1](https://github.com/Danderson123/Amira/compare/v0.9.0...v0.9.1) (2025-05-09)


### Documentation

* update E. coli panRG link ([d1199d8](https://github.com/Danderson123/Amira/commit/d1199d8a41429d57db7e6c9cf0e5f41f24ade83c))


### Bug Fixes

* limit graph building to 1 CPU to bypass slow runtimes across multiple CPUs, and update README with new E. coli panRG ([c93cabe](https://github.com/Danderson123/Amira/commit/c93cabe0c3c162503f0845121b2578518a838563))


### Styles

* pre-commit reformatting ([6abd6df](https://github.com/Danderson123/Amira/commit/6abd6dfb089c66e0632cdf00c59ea36af2d88f09))

## [0.9.0](https://github.com/Danderson123/Amira/compare/v0.8.0...v0.9.0) (2025-05-07)


### Continuous Integration

* do not test on python 3.9 ([ae1de18](https://github.com/Danderson123/Amira/commit/ae1de18cbc2f791417ff71ff63a21816157ec75a))
* do not test on python 3.9 ([d33391e](https://github.com/Danderson123/Amira/commit/d33391e1e13313dbbff7e26cc81bc20474664796))
* test on python 3.12 ([a6dc2c1](https://github.com/Danderson123/Amira/commit/a6dc2c149a5feada6d9fe77c57819fdac8951901))


### Miscellaneous Chores

* support S. pneumoniae and S. aureus ([31e0b0a](https://github.com/Danderson123/Amira/commit/31e0b0a81b2f394e189876b3c25f0d70e45a9aa4))
* support S. pneumoniae and S. aureus ([44eebb8](https://github.com/Danderson123/Amira/commit/44eebb84a5dee9d7c5659a51d4be3ad726224106))
* version bump ([93418c5](https://github.com/Danderson123/Amira/commit/93418c5fdd9cc4ccc60ff19fc69603ba8dc98d7c))
* version bump ([d9ffa24](https://github.com/Danderson123/Amira/commit/d9ffa2419fb93ecdff8052d46ddadbd57e44dcda))


### Documentation

* fix README ([03d6f40](https://github.com/Danderson123/Amira/commit/03d6f40818c8877867da1bb41b73a76e8d5dcc54))
* specify pandora version ([9fd68ab](https://github.com/Danderson123/Amira/commit/9fd68abd0b32c681b383a15380eba1d10fdd7e12))
* update download links in README ([20f278a](https://github.com/Danderson123/Amira/commit/20f278a7704abdc5043588ea7f99ad2510a878b9))
* update README ([46a2058](https://github.com/Danderson123/Amira/commit/46a2058c1a39bacdbc7bf9dcafa53dff587592bc))
* update README ([17c555b](https://github.com/Danderson123/Amira/commit/17c555b06b1046a0169072e29dfbf679a0c00bd0))
* update README ([27a6c37](https://github.com/Danderson123/Amira/commit/27a6c372b2d3424b3976aebde06314251da018b4))
* update README ([f630555](https://github.com/Danderson123/Amira/commit/f6305557365028424949e885f32aebe1fe6c6332))
* update README ([84af8be](https://github.com/Danderson123/Amira/commit/84af8be75cd2f8c8f5aa8bbe6b15eccdb85547a1))


### Features

* add SA and SP files ([fe041ef](https://github.com/Danderson123/Amira/commit/fe041ef5d372f082276e5d2fda0710cc037f98f3))
* increase subsample size ([0f0e353](https://github.com/Danderson123/Amira/commit/0f0e353927a30efd2c7562ef3999f141fe93b975))
* refactor graph correction to massively improve runtime ([70a8d66](https://github.com/Danderson123/Amira/commit/70a8d66e27cbf54aea58a3ed62b47fcdf325c1be))


### Bug Fixes

* allow reads with .fq extension ([38cf4db](https://github.com/Danderson123/Amira/commit/38cf4dba0e24771ca2895ad9131a5c7eede0b99b))
* allow reads with .fq extension ([8c3ac75](https://github.com/Danderson123/Amira/commit/8c3ac75a91fb5f3758c0cb0b96cb7912652ea1e1))
* allow underscores in read identifiers ([d2a3826](https://github.com/Danderson123/Amira/commit/d2a3826af38c328397fad3a3ed029f89b1e9af09))
* allow underscores in read identifiers ([d135d96](https://github.com/Danderson123/Amira/commit/d135d968ddd2a29a62d088c8e00d31e30a4328bd))
* attempt optimisation ([61cae8f](https://github.com/Danderson123/Amira/commit/61cae8f5887d6d81939e509648d27f7d38a7111c))
* bring up to speed with dev ([d50bb5f](https://github.com/Danderson123/Amira/commit/d50bb5f6117b2f892a53def33070aff578e18914))
* do not delete uncompressed fastq ([ea151b6](https://github.com/Danderson123/Amira/commit/ea151b6185e8752cb9c514498badc9a956699aee))
* do not write a modified fastq file ([59c013c](https://github.com/Danderson123/Amira/commit/59c013c28c57400d3d86732f57c4d21c50a85136))
* extend bubble popping path length ([6978c20](https://github.com/Danderson123/Amira/commit/6978c20ee283c54fb33671fdc491fc42e6851722))
* file path fix ([0a94737](https://github.com/Danderson123/Amira/commit/0a9473710c13996c3a4c2cec7f5c33df5cb5ba2c))
* fix bug in modifying read identifiers ([d309f47](https://github.com/Danderson123/Amira/commit/d309f473eebb0b92d759f2701b4de01e22d8a967))
* fix parsing of non-gzipped fastq ([4760b2c](https://github.com/Danderson123/Amira/commit/4760b2cbe0b23d402f92662237ee486dba2f45b0))
* fix unittests ([6610dc9](https://github.com/Danderson123/Amira/commit/6610dc962a5c192e26b14c32a66a7fec5699b0c8))
* further attempts to optimise bubble popping ([e43b0b5](https://github.com/Danderson123/Amira/commit/e43b0b55b787f726b59d89c3c1e8c39fb7d351ed))
* missing output dir ([b3d62d1](https://github.com/Danderson123/Amira/commit/b3d62d19d3a7d16da86d7c85c47fdc88ef407f24))
* missing output dir ([ffa58c0](https://github.com/Danderson123/Amira/commit/ffa58c0ad9b5d04517470461a3a4853dec76db97))
* new fastq path ([2de66e7](https://github.com/Danderson123/Amira/commit/2de66e7558a9ec733a9c9f06834dfd4ba84f60f6))
* new fastq path ([c3c4574](https://github.com/Danderson123/Amira/commit/c3c4574051539953ef8ab0e26e88db744e69a261))
* optimise path correction filtering ([be6dd67](https://github.com/Danderson123/Amira/commit/be6dd67e310bd378f9b5ddb834c8b89d9002ce49))
* optimise path filtering with suffix trees ([dd2a08a](https://github.com/Danderson123/Amira/commit/dd2a08aef097056115a7758ce857c0feb32d3bf6))
* optimise path filtering, leads to 10x speedup in some cases ([fb2e871](https://github.com/Danderson123/Amira/commit/fb2e8715753b5a5ccc8d8e04c51b6e0d5a496cd4))
* prevent unnecessary path comparisons ([4af472c](https://github.com/Danderson123/Amira/commit/4af472cc6568126516bd91fda3a390159b659421))
* progress reporting fix ([2781493](https://github.com/Danderson123/Amira/commit/2781493f267bc3306edd13687f7f5a732671ea44))
* reduce bubble popping path length slightly ([18caf38](https://github.com/Danderson123/Amira/commit/18caf38ca4c8210d964503868585ea75682a3472))
* reduce sourmash version ([2199d35](https://github.com/Danderson123/Amira/commit/2199d35b94da58abdaf963a6e5a1c782552e9c24))
* remove bounds from k-mer cutoff estimation ([128c8f2](https://github.com/Danderson123/Amira/commit/128c8f20552472ad3b6217583f2fc222d76c2f42))
* remove debug line ([d0b3c04](https://github.com/Danderson123/Amira/commit/d0b3c043a6bf6be73fbed16b096268b86258e4a0))
* remove debug message ([85c6713](https://github.com/Danderson123/Amira/commit/85c6713fce0eee352967b497713c22ea970e25f7))
* remove unnecessary progress report ([18c28d2](https://github.com/Danderson123/Amira/commit/18c28d2b499f31a25eac76cfdb24aa0c77134434))
* revert optimisation attempt ([22b5ec5](https://github.com/Danderson123/Amira/commit/22b5ec5af0164842595f46485f9eaa1c618c7f82))


### Styles

* pre-commit reformatting ([3db9704](https://github.com/Danderson123/Amira/commit/3db9704d6fd8e803c8ffae72a5e07ae129434ebd))
* pre-commit reformatting ([eab1299](https://github.com/Danderson123/Amira/commit/eab12999a4b7f1bba9639e4b4ba901a9739c35b6))
* pre-commit reformatting ([9283732](https://github.com/Danderson123/Amira/commit/928373253dd59940ce2caac8c8a9bf94277c7e43))
* pre-commit reformatting ([01cc43c](https://github.com/Danderson123/Amira/commit/01cc43cc170b78d8aeded6b7943e29bb6d4649fb))
* pre-commit reformatting ([78bc7fc](https://github.com/Danderson123/Amira/commit/78bc7fc5c6da1ae1c6c55770a7c26285a0ca462f))
* pre-commit reformatting ([cf0e509](https://github.com/Danderson123/Amira/commit/cf0e5098a22dd83110dda37e68da6a7789f0c8e2))
* pre-commit reformattng ([e1d7bf0](https://github.com/Danderson123/Amira/commit/e1d7bf01999e9c11ffdebfedd233d0336fbb4805))
* refactor main function ([e592795](https://github.com/Danderson123/Amira/commit/e59279512cb9933e4f22f9214128a59c1f71abb8))
* remove redundant reporting ([d1ad5a0](https://github.com/Danderson123/Amira/commit/d1ad5a0cd4843079aa3195b50aec1d74385f963c))
* rename cellular copy number column ([01adc54](https://github.com/Danderson123/Amira/commit/01adc5442c1220f606b2212b4b7a192b93c459de))

## [0.8.0](https://github.com/Danderson123/Amira/compare/v0.7.0...v0.8.0) (2025-03-25)


### Miscellaneous Chores

* version bump ([487c6b0](https://github.com/Danderson123/Amira/commit/487c6b0ce33cffb7062feb0d4a19a8cd65af5669))


### Documentation

* improve arguement documentation ([431dcab](https://github.com/Danderson123/Amira/commit/431dcabdeaeca49f0a741ea88d52a0b06eb17747))


### Features

* allow parsing of species files from command line ([a0dc317](https://github.com/Danderson123/Amira/commit/a0dc31700390a930fda3cc6d3e9cdcbecb45fa49))
* estimate copy numbers of paths and apply to covered AMR genes ([870f5e3](https://github.com/Danderson123/Amira/commit/870f5e361da297a124fe1c30ffcb771e237af1f4))
* use jellyfish to estimate cellular copy numbers ([5198fe2](https://github.com/Danderson123/Amira/commit/5198fe203729cdc342ca93a53f01e6ea624b631d))
* use kmer counts to estimate cellular copy number ([1ac9040](https://github.com/Danderson123/Amira/commit/1ac9040869091e09c9b36504ee7307bac72e92b7))


### Bug Fixes

* add jellyfish to singularity recipe ([00b24d4](https://github.com/Danderson123/Amira/commit/00b24d450bc435497122927273228e4700a66e86))
* add test files ([cee63ba](https://github.com/Danderson123/Amira/commit/cee63ba1eb7b017e4d3481da462a9f16a46616be))
* bug fix for non-redundant path filtering ([f5cfa38](https://github.com/Danderson123/Amira/commit/f5cfa38d643fb81db71b28d6ef890212cc8ffcb5))
* bug in gene positions when path does not need correcting ([a9a9652](https://github.com/Danderson123/Amira/commit/a9a9652b683fdb1ad3440b6102bf88f67b19e3f7))
* build container from main ([b65cc54](https://github.com/Danderson123/Amira/commit/b65cc5465b3f70d7760b0a262029cd34d12ffb63))
* correctly modify gene positions ([5925772](https://github.com/Danderson123/Amira/commit/5925772fddd8176e8b92171bf986924644486b2a))
* divide depth estimate by genomic copy number ([943577f](https://github.com/Danderson123/Amira/commit/943577f570c6e1e6c59d5f1280dbbd8733772e49))
* do not filter out whole paths ([adfa1eb](https://github.com/Danderson123/Amira/commit/adfa1eb2ee368b93ea54690beb94e265371b4809))
* don't specify jellyfish version in recipe ([818437d](https://github.com/Danderson123/Amira/commit/818437dead53c57d48194461bd47045dcc218112))
* estimate peak coverage from filtered histogram ([5c2ef2c](https://github.com/Danderson123/Amira/commit/5c2ef2c3864c89599395d6c351d277a102b49b4c))
* file extension bug ([429bfe6](https://github.com/Danderson123/Amira/commit/429bfe6844a1e52fdf2df94a3dd16b5c5b122d1e))
* filter alleles solely based on relative read depth ([f645a20](https://github.com/Danderson123/Amira/commit/f645a20dcf6c0ff6dc977d42fc4d644c6777104d))
* filter gene positions too ([a05b1ff](https://github.com/Danderson123/Amira/commit/a05b1ff207c0cdf08325fbaa874f953f444be880))
* filter genes based on relative frequency ([f559141](https://github.com/Danderson123/Amira/commit/f559141ba2262f0f9016b78672dfc6448b6cde5e))
* fix bug where path finding is blind to paths that start within an internal block ([0b0012d](https://github.com/Danderson123/Amira/commit/0b0012dbbcdd494aa8b9b1fae290a86d9911b157))
* fix bug where path finding is blind to paths that start within an internal block ([03f1020](https://github.com/Danderson123/Amira/commit/03f102036c145bd851092ad63639244b6844f53e))
* fix bug where path finding is blind to paths that start within an internal block ([f8dc458](https://github.com/Danderson123/Amira/commit/f8dc458ff27db994c249769e9c686c8869534937))
* gene coverage default val ([756870f](https://github.com/Danderson123/Amira/commit/756870f7d84b2a4180ae3d1ee3529ce6dd424c92))
* gunzip reads before giving to jellyfish ([c865567](https://github.com/Danderson123/Amira/commit/c86556770d58aec2f419f2f10e086ce9255aefd3))
* increase debugging information ([78e9408](https://github.com/Danderson123/Amira/commit/78e9408ad7e65fc5781942884d3002e9b499e9bf))
* increase debugging information ([87737e2](https://github.com/Danderson123/Amira/commit/87737e23099f1c580f14ffc1dd6c5c6d44c5a404))
* increase gene_min_coverage default ([ac38d91](https://github.com/Danderson123/Amira/commit/ac38d9196da1bc0fb826d70c72d96424f8d8b02a))
* increase k-mer size and filter by mean relative depth instead of copy number ([f879ce7](https://github.com/Danderson123/Amira/commit/f879ce774ca95e33ca83c4873ba1a10513d779c7))
* minor bug fix for read ids ([09b5d47](https://github.com/Danderson123/Amira/commit/09b5d4790d6b984aa21354fbdc6e13799b8f93db))
* ouput node coverage ([557f84a](https://github.com/Danderson123/Amira/commit/557f84a2fcba6f805827b0976636e8d026ce4a6a))
* reduce kmer size ([f78ab56](https://github.com/Danderson123/Amira/commit/f78ab561b73f470a7206f8c7c939dbe436aa865b))
* remove print statement ([9fc6349](https://github.com/Danderson123/Amira/commit/9fc63498e56ece9530cd64412aee9d28e276735c))
* remove redundant tests ([b5da866](https://github.com/Danderson123/Amira/commit/b5da866d1d5dcecd714a3c4d0b759c15c6d5c73c))
* remove redundant variables ([08230cd](https://github.com/Danderson123/Amira/commit/08230cd006e3edafbe68ece2837443cc9a8f5028))
* removes contamination filtering and replaces with filtering by depth ([428d01a](https://github.com/Danderson123/Amira/commit/428d01a9f8479f52adbd2682e713ef14b70b4cc5))
* resolve conflicts ([ac6c46f](https://github.com/Danderson123/Amira/commit/ac6c46f68ecd4a7bf0f53ec45ef3bc29737dcd5f))
* revert back to short reads ([f39e422](https://github.com/Danderson123/Amira/commit/f39e4225fb649bd70d6752bbbebbee113cf8f9bf))
* skip relative depth based filtering if read coverage is low ([eb83fad](https://github.com/Danderson123/Amira/commit/eb83fad1b0dfc9241f56261bcfdcbe80e7335a10))
* skip supplementary reads ([6115e40](https://github.com/Danderson123/Amira/commit/6115e40a6668e83d711e3e6d90b522e32b55308a))
* stop pandora from filtering based on coverage ([f77d228](https://github.com/Danderson123/Amira/commit/f77d228ccd6758f66807eb4d8fd6bd89daa1741f))
* unittest fix ([6f006a6](https://github.com/Danderson123/Amira/commit/6f006a6ce23a2998dbe54d7fb5c3cfbecaa52680))
* unittest fix ([05bd1ec](https://github.com/Danderson123/Amira/commit/05bd1ec91f5a37863d4b7e5b439806cd0d0c2193))
* use estimated read depth when kmer depth estimate is too high ([d4a411f](https://github.com/Danderson123/Amira/commit/d4a411f12959b39a03265143c5de75d7da623a23))
* use initial reads to supplement AMR genes ([e128d57](https://github.com/Danderson123/Amira/commit/e128d579a9b866a4f1530c283363ab462684b7f6))
* use kmer depth to get copy number estimate ([bfeddba](https://github.com/Danderson123/Amira/commit/bfeddba8db1faab18ecbb888be66ec3cfd4da586))
* use main branch for Singularity container ([899311f](https://github.com/Danderson123/Amira/commit/899311f568dac64f3fbc72f61500e5b79a7f244d))
* use scipy find peaks to estimate mean k-mer depth ([a85faf4](https://github.com/Danderson123/Amira/commit/a85faf42eb05013703ac595f213b3af2198d1496))


### Styles

* more consistent stderr messages ([b4a72b4](https://github.com/Danderson123/Amira/commit/b4a72b49806e7b0ec62224cd63c528b39131de72))
* pre-commit reformatting ([2ca3bc1](https://github.com/Danderson123/Amira/commit/2ca3bc1ec070202f7386bc5426eef03cd94df815))
* pre-commit reformatting ([697e0e3](https://github.com/Danderson123/Amira/commit/697e0e344b1f95dd27ada1bb743f8c2c0d15e226))
* pre-commit reformatting ([4fd4030](https://github.com/Danderson123/Amira/commit/4fd403070a08306fd89dd0c536ee11e4ddc941d2))
* pre-commit reformatting ([01de505](https://github.com/Danderson123/Amira/commit/01de505f909aa4d708090c006a074a8eda053110))
* pre-commit reformatting ([86bab63](https://github.com/Danderson123/Amira/commit/86bab63d7ff6c79b5ea14215d350751ae632b4c6))
* pre-commit reformatting ([e06c36d](https://github.com/Danderson123/Amira/commit/e06c36d7f3bd1b28d4df7706ac66ba181235e82f))
* pre-commit reformatting ([8084235](https://github.com/Danderson123/Amira/commit/8084235471d7f380493f5971a2a5bcf3eee43afd))
* pre-commit reformatting ([7270d99](https://github.com/Danderson123/Amira/commit/7270d99c9a30a4e7775a836fe78248ad98ad8be7))
* pre-commit reformatting ([7182e75](https://github.com/Danderson123/Amira/commit/7182e757231ee7b96dd99475e1029ef81a1d2ebc))

## [0.7.0](https://github.com/Danderson123/Amira/compare/v0.6.5...v0.7.0) (2025-03-04)


### Miscellaneous Chores

* version bump ([2ac3804](https://github.com/Danderson123/Amira/commit/2ac3804f88e261d95770c7075122511fde9ff362))


### Documentation

* add download link to E. faecium panRG ([65e77e5](https://github.com/Danderson123/Amira/commit/65e77e57cab4b839056e9b892f8e4463f13ed90e))


### Features

* do not filter suspected contaminants by default ([86c71fe](https://github.com/Danderson123/Amira/commit/86c71fe952ddd79d40d134c08f5324b65d5a58eb))
* support Enterococcus faecium ([3bfd3a3](https://github.com/Danderson123/Amira/commit/3bfd3a3c19e9b6f32d5c98b06bbf10716baa2a52))

## [0.6.5](https://github.com/Danderson123/Amira/compare/v0.6.4...v0.6.5) (2025-02-26)


### Miscellaneous Chores

* version bump ([1792fa4](https://github.com/Danderson123/Amira/commit/1792fa4306a3a1c8f4b7566822d1ae62f56265ac))


### Bug Fixes

* reduce pandora gene filtering coverage ([ecc0822](https://github.com/Danderson123/Amira/commit/ecc082243a48c881b68f840c20119834d5ee759a))
* remove special characters from read identifiers ([1088e2c](https://github.com/Danderson123/Amira/commit/1088e2c2602eb07678d6f28369f23fb77be5522e))
* write fastq with updated read identifiers to output directory ([dc5c197](https://github.com/Danderson123/Amira/commit/dc5c1971475be0516083b6be77681ebb66c2a785))


### Styles

* pre-commit reformatting ([07fdd57](https://github.com/Danderson123/Amira/commit/07fdd575e19e2bd13be1c78ce6cbcbf38504f806))

## [0.6.4](https://github.com/Danderson123/Amira/compare/v0.6.3...v0.6.4) (2025-02-13)


### Miscellaneous Chores

* version bump ([aef08a3](https://github.com/Danderson123/Amira/commit/aef08a31d06d38e843dac8a659c8275fbf8b15f7))


### Bug Fixes

* increase polishing iterations ([a409a12](https://github.com/Danderson123/Amira/commit/a409a1233db72a9f57927ed2c7bd5a0495913fd3))
* increase window size to prevent a racon bug from incorrectly adding bases to alleles ([08bd3b2](https://github.com/Danderson123/Amira/commit/08bd3b29c03441d7695a61342ae83b3c16dd58c2))
* tweak decision process for alleles ([c0652b4](https://github.com/Danderson123/Amira/commit/c0652b4f299969677ecefb3b912c76bb6e30b23d))


### Styles

* pre-commit reformatting ([93c19ef](https://github.com/Danderson123/Amira/commit/93c19ef951c2947fcf34c6ff731d61d1772ce280))

## [0.6.3](https://github.com/Danderson123/Amira/compare/v0.6.2...v0.6.3) (2025-02-12)


### Miscellaneous Chores

* version bump ([9165b0a](https://github.com/Danderson123/Amira/commit/9165b0a599a320dea44bfb81ecd7dc92ce7540d8))


### Bug Fixes

* bug fix for coverage filtering ([6b1890f](https://github.com/Danderson123/Amira/commit/6b1890f3a5739e1014d477ba224fe82cebf7d418))
* bug fix for read sampling ([7a05734](https://github.com/Danderson123/Amira/commit/7a0573434c1672da3afd720dd3ded162bdadfba1))
* bug fix in path finding when duplicated nodes occur next to eachother ([48eedad](https://github.com/Danderson123/Amira/commit/48eedad3552f5ba2c679d8dee6208ef07ac94e40))
* build container from main branch ([87217e6](https://github.com/Danderson123/Amira/commit/87217e60d3e520ced79ab1644b20ec2bd85cc76f))
* compensate for unmappped reads with pandora ([8f392a8](https://github.com/Danderson123/Amira/commit/8f392a8be97f8fb996ecbbbb7e79a49b159827bf))
* prioritise alignment coverage when finding closest allele ([06aca6c](https://github.com/Danderson123/Amira/commit/06aca6c1f0be3d8a096c575402a091055cf98f91))
* remove redundant file ([90a7aee](https://github.com/Danderson123/Amira/commit/90a7aee4277b24f2e3e4ac371894d41694f7a786))
* run pandora on full read set before sampling ([27051a3](https://github.com/Danderson123/Amira/commit/27051a355395083b86d2bc6d757bdb2e56f38101))
* run poetry lock before installation ([83fbfd1](https://github.com/Danderson123/Amira/commit/83fbfd15dd743ef4b21ead7000f5b1d965e0271f))
* seed pandora to ensure runs are reproducible ([2f00f98](https://github.com/Danderson123/Amira/commit/2f00f987d00b428b1b198164b66926a5a4f7aab7))
* sort dictionary of read annotations before sample ([275f57f](https://github.com/Danderson123/Amira/commit/275f57f6d0fa3234bc232c88d7f8f6f917b5eb42))
* sort dictionary of read annotations before sampling ([2acde47](https://github.com/Danderson123/Amira/commit/2acde470a923299a7f947c5a2cf6fb9b70fc3b42))
* update __main__.py ([f98e3ea](https://github.com/Danderson123/Amira/commit/f98e3ea0fc5375bc37e3a27da9ec7a401cd96727))


### Styles

* pre-commit reformatting ([11028f2](https://github.com/Danderson123/Amira/commit/11028f2c97f56914775eb18f4897472f3fc0d663))
* pre-commit reformatting ([1f92286](https://github.com/Danderson123/Amira/commit/1f92286243cf3a6a35d5d0f885957f3d0ef61e80))
* pre-commit reformatting ([1a792d8](https://github.com/Danderson123/Amira/commit/1a792d861ac155a363db14cbd0f4263549403ff6))

## [0.6.2](https://github.com/Danderson123/Amira/compare/v0.6.1...v0.6.2) (2025-02-09)


### Miscellaneous Chores

* version bump ([8f97ebc](https://github.com/Danderson123/Amira/commit/8f97ebc60084345723799dce531523780d41ff0c))
* version bump ([ab17ff1](https://github.com/Danderson123/Amira/commit/ab17ff13b5b8b1d0788a5d747273ea9169383d39))

## [0.6.1](https://github.com/Danderson123/Amira/compare/v0.6.0...v0.6.1) (2025-02-09)


### Bug Fixes

* bug fix when alleles filtered ([fe5ccd4](https://github.com/Danderson123/Amira/commit/fe5ccd4718ba6a5600d4307fd8ff36b063a4cf2d))

## [0.6.0](https://github.com/Danderson123/Amira/compare/v0.5.0...v0.6.0) (2025-02-09)


### Build System

* add singularity recipe and subsample reads before calling pandora ([a1d16a1](https://github.com/Danderson123/Amira/commit/a1d16a1c4bfb44fce49786f014327bdef9bb9edc))


### Miscellaneous Chores

* remove singularity build from release CI ([33550d2](https://github.com/Danderson123/Amira/commit/33550d2273c2f16ec53ab4a4535ae4c90c73ddb6))
* remove singularity build workflow ([e173b21](https://github.com/Danderson123/Amira/commit/e173b21b45cc436d826e6611d54a099570867729))
* singlarity build yaml ([c6217d9](https://github.com/Danderson123/Amira/commit/c6217d94e2f7ec707c7c47d12d49be31d80afa63))
* update version ([4accc24](https://github.com/Danderson123/Amira/commit/4accc2471de68d5a2e71c9e203dec73a817f548e))
* update version ([4ef0302](https://github.com/Danderson123/Amira/commit/4ef0302168816ba5d1edd28e61816844d00dca66))


### Documentation

* add download link for E. coli panRG ([31300d0](https://github.com/Danderson123/Amira/commit/31300d0f64c70473fe3245e58b757e74265bf466))
* add download link for K. pneumoniae panRG ([f8e4bd3](https://github.com/Danderson123/Amira/commit/f8e4bd3ed24a72cbc7dd0aa233eea95eaa6ec82a))
* add instructions to build singularity container ([4c35f60](https://github.com/Danderson123/Amira/commit/4c35f60be44436f575cfc36dbd9795f672cf891d))
* modify for clarity ([ff96ba5](https://github.com/Danderson123/Amira/commit/ff96ba5d3d75c31b732eed89e9c0f2182b938c2e))
* update README ([2349309](https://github.com/Danderson123/Amira/commit/2349309956a8b848f14401e0cb5be2fec237c7c6))


### Features

* add reference files for K.pneuomoniae ([29146a1](https://github.com/Danderson123/Amira/commit/29146a1abe8df45a920d48e4f0d28bd80d70f4f3))
* add singularity container build to CI ([c7d9b11](https://github.com/Danderson123/Amira/commit/c7d9b115f0e93032a7ea8f562970a59d9c92cd7d))
* filter contaminants by default ([9ee5fa1](https://github.com/Danderson123/Amira/commit/9ee5fa1c21b59b4394160884abcacb794af1db7f))
* integrate pandora into amira and specify species from command line ([e3fb93c](https://github.com/Danderson123/Amira/commit/e3fb93c7e49518afc2aef8cb976a1e4fc4e1d0bc))
* make some arguements default and allow specification of AMR gene coverage and identity. ([5bae165](https://github.com/Danderson123/Amira/commit/5bae16540752a270f3d5533853f89394e6c30522))
* randomly subsample reads by default to improve runtime ([bc28cab](https://github.com/Danderson123/Amira/commit/bc28cab0c499deac7cee23d3533d09a84d2a1350))
* subsample reads before parsing pandora sam ([12e3d0b](https://github.com/Danderson123/Amira/commit/12e3d0bed2fcf952599c314983e21464da98b5e9))


### Bug Fixes

* add samtools path variable ([c5f5e16](https://github.com/Danderson123/Amira/commit/c5f5e167e0ba575ab6c35bca734ff2b51d4eb78d))
* correctly access species files from module ([b232179](https://github.com/Danderson123/Amira/commit/b232179113d3c910f4680eaeaacec7b211a372ec))
* fix singularity recipe ([961c23c](https://github.com/Danderson123/Amira/commit/961c23c9ffe2e701fbb6c39ab7a6c722fd09efb3))
* fix singularity recipe ([7bfaca1](https://github.com/Danderson123/Amira/commit/7bfaca19a374cb8899880ea3365563f6e923f00c))
* fix singularity recipe ([2622171](https://github.com/Danderson123/Amira/commit/26221718c8679ca40a59424fb75ce44188c3c0c1))
* fix singularity recipe ([f32fd59](https://github.com/Danderson123/Amira/commit/f32fd59c45c0d2f4cdcae6127877fdc1f5b70963))


### Styles

* pre-commit reformatting ([094bfdc](https://github.com/Danderson123/Amira/commit/094bfdcb147485a26335a4a9537b71c0ad437a3b))
* pre-commit reformatting ([63b54f3](https://github.com/Danderson123/Amira/commit/63b54f39319fc2899d987a11a2b3e7950b470407))
* pre-commit reformatting ([23bd536](https://github.com/Danderson123/Amira/commit/23bd536a7ff65af7db54d16a992fb9b2fca2c6f9))
* pre-commit reformatting ([884228d](https://github.com/Danderson123/Amira/commit/884228d284b027bba2510eeabfd6f9d8862cf730))


### Tests

* update unittests ([e544a1f](https://github.com/Danderson123/Amira/commit/e544a1fdfef868c802d7b2bee54b46618015579b))

## [0.5.0](https://github.com/Danderson123/Amira/compare/v0.4.1...v0.5.0) (2025-01-24)


### Continuous Integration

* pre-commit reformatting ([4834ee0](https://github.com/Danderson123/Amira/commit/4834ee0e87c1ae0a63f4e727f82e794236de044d))


### Miscellaneous Chores

* fix version ([6967db9](https://github.com/Danderson123/Amira/commit/6967db9c4b218f05869a6a8cec7e961a28f8ad3e))
* version bump ([58cdd6d](https://github.com/Danderson123/Amira/commit/58cdd6dd22d8b8c95f9fafb912504b6c1d34b076))


### Documentation

* update README ([08458a6](https://github.com/Danderson123/Amira/commit/08458a6dfd281b4191bfde8366a775d27239cc40))


### Features

* detect AMR-determining SNPs in promoters if specified ([582ea7a](https://github.com/Danderson123/Amira/commit/582ea7ae5ab3b9b037cd00632c39859a4f41a76c))
* estimate cellular copy numbers by calculating coverage across the longest read for each AMR gene, then normalise by mean read depth across core genes ([68566de](https://github.com/Danderson123/Amira/commit/68566def6ddbe86a4f8a1ce36cb527b313666868))
* genotype novel SNPs in promoters ([ba8e654](https://github.com/Danderson123/Amira/commit/ba8e654fbbc20bf1b0ec497509ad9bed167f43b2))


### Bug Fixes

* check for float ([0eef5b2](https://github.com/Danderson123/Amira/commit/0eef5b2c2a244f1e7cddda4ad70c0cd121d46e6d))
* detect presence of causative SNPs without full allele matches] ([7512592](https://github.com/Danderson123/Amira/commit/75125926c088a0975a761388e704ef056f443fd4))
* estimate mean depth only from core genes ([49006e4](https://github.com/Danderson123/Amira/commit/49006e4fc6bec702d85b979f009808717e422694))
* fastq file path fix ([52a6fbb](https://github.com/Danderson123/Amira/commit/52a6fbb06707012015637679a09618949594abc4))
* fix adding promoter results to existing df ([2a55d4f](https://github.com/Danderson123/Amira/commit/2a55d4f42605adb165af3d94055976d2e0477cb0))
* fix adding rpomoter results to existing df ([8fa1c92](https://github.com/Danderson123/Amira/commit/8fa1c9243bf522f1d67f2a4ed900336c7fa0b5e9))
* only add promoters to results if non ref ([6afa85b](https://github.com/Danderson123/Amira/commit/6afa85b197c3ecbcddc230e7a81d70c0926f679b))
* only add promoters to results if non ref ([1c30eec](https://github.com/Danderson123/Amira/commit/1c30eec5d581adbc4d8337e123846d680b812b8e))
* remove depreciate append to df ([b0f63e9](https://github.com/Danderson123/Amira/commit/b0f63e91f587c34a72532bda2bde274b7971e5a2))
* remove underscores in read names and remove redundany graph outputs ([c28172c](https://github.com/Danderson123/Amira/commit/c28172c28ba21f347e801e4206775a134a42680c))
* tweak mpileup params ([f23bfa2](https://github.com/Danderson123/Amira/commit/f23bfa2c70d2d7b1e9b4612e201e9e79ad27fbd5))
* update graph_operations ([921953b](https://github.com/Danderson123/Amira/commit/921953b38595f26e2d2828d15c2cd69c5fc786bd))
* use min coverage of node in path to estimate copy number ([d5c1291](https://github.com/Danderson123/Amira/commit/d5c1291f1bdc5168f87d5da33c22e8d12c7396c8))
* use samtools mpileup to estimate cellular copy number to include secondary alignments ([440096e](https://github.com/Danderson123/Amira/commit/440096e9437d209dfffd418f081eee60b6e75ea6))


### Styles

* flake8 formatting ([cde598d](https://github.com/Danderson123/Amira/commit/cde598d2ce263f892bec5f322469cb662f3b4b1f))


### Tests

* add test files ([d698e0a](https://github.com/Danderson123/Amira/commit/d698e0a1fa815a104d08ae8a40b0557c342b8675))

## [0.4.1](https://github.com/Danderson123/Amira/compare/v0.4.0...v0.4.1) (2024-12-16)


### Bug Fixes

* use node-space path to estimate copy number in context ([1f58801](https://github.com/Danderson123/Amira/commit/1f58801ecccbec8347102009bc09e3a2df8a7e67))

## [0.4.0](https://github.com/Danderson123/Amira/compare/v0.3.1...v0.4.0) (2024-12-15)


### Miscellaneous Chores

* version bump ([99a2389](https://github.com/Danderson123/Amira/commit/99a23893db7efa69fd713b276b592ca4080a9622))


### Features

* use gene-space paths to assign reads ([b30ede6](https://github.com/Danderson123/Amira/commit/b30ede6a2f4faa1ca3306c2172514faf23f479c0))


### Bug Fixes

* add singleton paths ([d7fe631](https://github.com/Danderson123/Amira/commit/d7fe6318f2cef55b2032803a4dcb35a90a248c48))
* apply new condition to gene-space subpath generation ([2fa0755](https://github.com/Danderson123/Amira/commit/2fa075515f3df2d8f54cfc0736f5b50b30bc6663))
* prevent checking of redundant path sublists ([f79924e](https://github.com/Danderson123/Amira/commit/f79924ea16f22b4c8f1c0b3012679e81e1946dc9))
* prevent combinatorial increase of gene-space sublists of paths. Multi-process with mp.pool ([0afcba5](https://github.com/Danderson123/Amira/commit/0afcba5a521f1adb1b3cabf1ee1af59240586b8a))
* prevent removal of small, high coverage components ([870a07e](https://github.com/Danderson123/Amira/commit/870a07e1ccadeaffd0c8e219fe03d34df5855aeb))


### Tests

* update unittest ([e211da2](https://github.com/Danderson123/Amira/commit/e211da2649bac7431d57c0ff3b70551a29f9fae1))
* update unittests to reflect changes ([f491c92](https://github.com/Danderson123/Amira/commit/f491c929fb42a2aeaec70aa0a953108690c2ec9d))

## [0.3.1](https://github.com/Danderson123/Amira/compare/v0.3.0...v0.3.1) (2024-11-18)


### Build System

* **deps:** bump tqdm from 4.64.1 to 4.66.3 ([b4d9b4f](https://github.com/Danderson123/Amira/commit/b4d9b4fb9814a0bb343fba395451cf95a4c1a9fe))
* **deps:** bump tqdm from 4.64.1 to 4.66.3 ([ee877f1](https://github.com/Danderson123/Amira/commit/ee877f1d1c4ef0732af88bea88b4d7e358787416))

## [0.3.0](https://github.com/Danderson123/Amira/compare/v0.2.0...v0.3.0) (2024-11-18)


### Miscellaneous Chores

* version bump ([04b1bef](https://github.com/Danderson123/Amira/commit/04b1beffc4618dc42c33b103a56dbb3c70242185))


### Documentation

* add citation placeholder ([0e1168e](https://github.com/Danderson123/Amira/commit/0e1168e6347d47a9386f32b79d2c538afebb44d6))


### Features

* use suffix trees to obtain all paths that span from the start of a block of AMR genes to the end ([db6eb3c](https://github.com/Danderson123/Amira/commit/db6eb3c717227653485084b7c8a943261a0be31c))


### Bug Fixes

* supplement full paths with singleton nodes ([7f6e028](https://github.com/Danderson123/Amira/commit/7f6e028ce44ab3b1a3c373480e0c45d36f3b90af))


### Code Refactoring

* import modularity of graph operations ([a2ebde5](https://github.com/Danderson123/Amira/commit/a2ebde5ddf53e7dd7b56d3f928548d58c180f047))
* improve modularity of read operations ([425370a](https://github.com/Danderson123/Amira/commit/425370aaa52b344a8ee8ae87f307c50451ad89a7))
* improve modularity of result operations ([747276f](https://github.com/Danderson123/Amira/commit/747276f47ae5d117e69ee3a3a6faee7544ec18f4))


### Styles

* pre-commit reformatting ([580064a](https://github.com/Danderson123/Amira/commit/580064a5456c3dd1b2056744944c28e0d14f4cc4))


### Tests

* add missing test files ([aa552a1](https://github.com/Danderson123/Amira/commit/aa552a18dfa760c7a7fc6f34bcd2d6064f0e40c6))

## [0.2.0](https://github.com/Danderson123/Amira/compare/v0.1.0...v0.2.0) (2024-11-06)


### Build System

* version bump ([eb97b98](https://github.com/Danderson123/Amira/commit/eb97b98307fd2f24abf33c3c4666860e95504d6d))


### Documentation

* update README ([162fb3b](https://github.com/Danderson123/Amira/commit/162fb3bcdc191c6059052dfc24e45ef3ec2cb6b8))
* update README ([0d87fc7](https://github.com/Danderson123/Amira/commit/0d87fc7002e15ce0b60761071fe7d0f7df720d4f))
* update README ([22be16a](https://github.com/Danderson123/Amira/commit/22be16a22807b92ad79aca166d569af0d08da1a4))


### Features

* display amira version ([bf241fb](https://github.com/Danderson123/Amira/commit/bf241fbdece7f0595b9cda5820a4765e475dee0d))


### Bug Fixes

* add version to init and pyproject ([eeecf4f](https://github.com/Danderson123/Amira/commit/eeecf4fb336b460128c662e6e3395b9d6f86fc53))
* bypass subsampling from dictionary ([ac7604e](https://github.com/Danderson123/Amira/commit/ac7604eb2687b2225c6f291bc58abad3c1b969c5))


### Styles

* pre-commit recommendations ([6b808e1](https://github.com/Danderson123/Amira/commit/6b808e10dc9e9c090d1e0ce3ef10eb40e170548d))

## 0.1.0 (2024-11-06)


### Continuous Integration

* add .isort.cfg ([47a69b5](https://github.com/Danderson123/Amira/commit/47a69b572301749e708be8248fb81a3039b45be7))
* add pre-commit to dependencies ([608333e](https://github.com/Danderson123/Amira/commit/608333ed8e7ae9f95ad0563e650edc4014135259))
* adds ci.yaml to run makefile ([ad91b87](https://github.com/Danderson123/Amira/commit/ad91b875a7d6974562349d9c602e2ab837f88a60))
* adds poetry.lock, .pre-commit-config.yaml and Makefile ([73810a3](https://github.com/Danderson123/Amira/commit/73810a32796660135766b8bf6080ff5b4e23d938))
* change release-please repo name ([f6126ef](https://github.com/Danderson123/Amira/commit/f6126efd5a702a75a58f33ce718f322ac4fdb425))
* comment out redundant tests ([f0ad490](https://github.com/Danderson123/Amira/commit/f0ad4900830a7f06827ac37123cb2b9eb3868371))
* comment out redundant tests ([1b41b3c](https://github.com/Danderson123/Amira/commit/1b41b3c90eec16f8060e1c8e72259acbc47f9425))
* create release-pypi.yaml ([3b9ae5c](https://github.com/Danderson123/Amira/commit/3b9ae5c7491a55018253aed253e3479c4b60dedb))
* Creates release-please.yaml ([53abf05](https://github.com/Danderson123/Amira/commit/53abf057f04a2949ee8f3f7879b5495fba5ec389))
* install dependencies before test ([7fc194a](https://github.com/Danderson123/Amira/commit/7fc194a8fac8e265e158e5bd3dfa529935c9db42))
* install dependencies before test ([2d13af7](https://github.com/Danderson123/Amira/commit/2d13af725b934f4f1b24d72b053d3a2e2be6fc1d))
* install dependencies in actions ([be2a163](https://github.com/Danderson123/Amira/commit/be2a163ddc85af89519503ad09c0cd19f3487809))
* install pytest before running tests ([79e8011](https://github.com/Danderson123/Amira/commit/79e8011ed5057e80e98e746fe7fee8a3ee24ecdc))
* Merge branch 'CI_integration' of https://github.com/Danderson123/amira_prototype into CI_integration ([44506b2](https://github.com/Danderson123/Amira/commit/44506b2a31e1ae77f3307aa67dd55fcc4888c01d))
* remove python3.8 support ([01d6a0e](https://github.com/Danderson123/Amira/commit/01d6a0e02701c5db6dc3338d5380bddc1f99d807))
* remove unused test script ([c8947f5](https://github.com/Danderson123/Amira/commit/c8947f5169396737190dab1ad1668697ac0ac56e))
* remove Windows support ([466cbcf](https://github.com/Danderson123/Amira/commit/466cbcf14a164bd5526176bc4f965bb8541c8b8b))
* removes redundant .github/workflows/run_tests.yml ([869fd5d](https://github.com/Danderson123/Amira/commit/869fd5d31552f01a2eeb56b89c6601519f72705e))
* set up python installation ([278a4a5](https://github.com/Danderson123/Amira/commit/278a4a511e9f8cecafb3052f686eefb2a3b89f9a))
* update __main__.py ([2cdeebb](https://github.com/Danderson123/Amira/commit/2cdeebb2a6affaab75072586cc05627176738f4a))
* update ci.yaml ([028b3c3](https://github.com/Danderson123/Amira/commit/028b3c3ad879cc0ec8a470c845834a57b39c324e))
* update ci.yaml ([bda10f6](https://github.com/Danderson123/Amira/commit/bda10f62f21abe4441c77aaab7a56a6672877d76))
* update ci.yaml ([f74ee4f](https://github.com/Danderson123/Amira/commit/f74ee4f5adc71ceecea01c7c13d589975361e1c3))
* update ci.yaml ([95a4c46](https://github.com/Danderson123/Amira/commit/95a4c466c2bcfb06120bc1fd3167a5d2a9c47bca))
* update poetry.locl ([1ab43b4](https://github.com/Danderson123/Amira/commit/1ab43b4c27f271e99ec2b264bc3c93368d81e252))
* update pre-commit tools using poetry ([40e7adc](https://github.com/Danderson123/Amira/commit/40e7adc6b662ae73e95dde37ad73ea4b94ac7e04))
* update release-please.yaml ([30e89e4](https://github.com/Danderson123/Amira/commit/30e89e42cb93754fd7864ac86b074aec8f8aec23))
* use pip to install pypi dependencies ([392864f](https://github.com/Danderson123/Amira/commit/392864fef9c9d40363e12f798c4459016080a26f))


### Documentation

* update README ([82e282b](https://github.com/Danderson123/Amira/commit/82e282be93903114a3d0957dbef1e8096f130471))
* update README ([00d8263](https://github.com/Danderson123/Amira/commit/00d826392cff7f4a60d3704f98edae0d4dc95aab))
* update README ([056c60c](https://github.com/Danderson123/Amira/commit/056c60c589e155e47e2d073985e0a547d49a2ff0))
* update README ([a178cab](https://github.com/Danderson123/Amira/commit/a178cab32d09d3babca819ff43ef8e3419d72b1f))
* update README ([2e427d6](https://github.com/Danderson123/Amira/commit/2e427d69972f43fcd0950f48f2aee434a4b41a14))
* update README and add ref files ([a6f85ae](https://github.com/Danderson123/Amira/commit/a6f85ae5446f1b488c99389f9b9c2a1e6091cb74))


### Features

* add closest reference allele header to read clustering keys for evaluation ([5967028](https://github.com/Danderson123/Amira/commit/5967028b17ab0c605c12345387d8659b259ca0db))
* add dynamic determination of node threshold ([307744c](https://github.com/Danderson123/Amira/commit/307744c36b31ca81635367b5e4c13c48e13dabf1))
* add partial untittests for new functions ([4eb534b](https://github.com/Danderson123/Amira/commit/4eb534b09935dd90d50ebc7a836e0e65ecf6f5bb))
* add source direction to gml ([63d3b38](https://github.com/Danderson123/Amira/commit/63d3b38d21e436a4bced929b05edea77c9bc72af))
* applies a fixed path threshold for correction  then uses minimizer containment of paths to correct high coverage paths in the final round. ([824bec0](https://github.com/Danderson123/Amira/commit/824bec0a59baafaddeaa19fadda15be6508027f0))
* correct bubbles until a maximum of 10 iterations ([7d4c87e](https://github.com/Danderson123/Amira/commit/7d4c87ebd320da91f0ceeb6806eef7cf834d8434))
* dynamic determination of post-bubble popping node filtering thresholds ([e84bc97](https://github.com/Danderson123/Amira/commit/e84bc97eb49b827e655652b4f97dd66ba6ece8ba))
* dynamically choose a value for k ([dfee11a](https://github.com/Danderson123/Amira/commit/dfee11a4f1c45b84c8af9296395678122661ecf6))
* filters alleles that are not &gt;=90% similar to a reference allele and polishes reference alleles to obtain Amira allele ([9d3e0e4](https://github.com/Danderson123/Amira/commit/9d3e0e41cc2665588d0cdffa8572bb3147d52ff1))
* initial attempts to correct all bubbles at once ([cd5fbae](https://github.com/Danderson123/Amira/commit/cd5fbae2cd4ff25ce7d6e495044c74e81c88b865))
* initial attempts to correct all bubbles at once ([15c40cd](https://github.com/Danderson123/Amira/commit/15c40cd85ce54503e9793006509783ae0adca872))
* kind of working correction method ([e79d6cf](https://github.com/Danderson123/Amira/commit/e79d6cf076c0d8a0dbea66d837d83ec115cfb738))
* mulit-processed DFS of paths and adds minimizer comparison of paths for accurate correction ([9718b80](https://github.com/Danderson123/Amira/commit/9718b801974e48f4e6577212ad3460fa8be47fe5))
* multiprocessed graph building using joblib ([c6af7ac](https://github.com/Danderson123/Amira/commit/c6af7ac427b2ba6835cc23c56ff47b1194a04041))
* new approach to cluster reads based on paths through the graph ([0c36248](https://github.com/Danderson123/Amira/commit/0c36248a12c3e61f69a91ea165925fbae9a3d7c9))
* output approximate copy numbers and reference allele depths ([fbe147b](https://github.com/Danderson123/Amira/commit/fbe147b4149712c74150f6d79a798766a9ff2fcf))
* output tsv of amira results ([cbe3c96](https://github.com/Danderson123/Amira/commit/cbe3c960d636cf122012eb7d855e382e5a8a628e))
* output txt file of depths across reference alleles ([5128c5c](https://github.com/Danderson123/Amira/commit/5128c5cd2aeb8cf2edd3bf87c4d23ee5d8c287ad))
* recover reads from paths that have been filtered out ([ef58223](https://github.com/Danderson123/Amira/commit/ef58223afcdfeab171e79fac94b89d66c3ab76cd))
* report all equally close variants ([d373b7e](https://github.com/Danderson123/Amira/commit/d373b7e7af356940881d8b28427770d3d7e0f7ae))
* take fasta reference file as input and polish closest AMR gene to get nucleotide sequence ([3255f28](https://github.com/Danderson123/Amira/commit/3255f28fea6b9295eaa8245d0a24b698cf9a2696))
* ten iterations of racon polishing ([a81ddce](https://github.com/Danderson123/Amira/commit/a81ddce14cfafbc2e4e1e53ead0b85be537c728a))
* tracks the position of each gene through each correction step to allow extraction of its sequence. Also removes redundant code and tests and improves code modularity. ([5a21e91](https://github.com/Danderson123/Amira/commit/5a21e91b1c844730f3aece97d54d1b37cb26bfee))
* use allele sequences to polish pandora consensus instead of entire read ([cd07418](https://github.com/Danderson123/Amira/commit/cd0741859598d950c65bea8d50b3f44547ed0dad))
* use as many adjacent nodes as possible to resolve complex AMR gene paths ([ca677c5](https://github.com/Danderson123/Amira/commit/ca677c5144c0223c80bbf51e575c4d529ec46f94))
* use overall mean node coverage for specific k to estimate copy number ([ece0037](https://github.com/Danderson123/Amira/commit/ece0037f8909532002677c832a76ea54e5426d56))


### Bug Fixes

* accounts for racon runtime errors, prevents AMR gene bubbles being popped ([d0be7f4](https://github.com/Danderson123/Amira/commit/d0be7f4bd1158c75b2a795bccc228283d642cca1))
* add dependencies to pyproject.toml ([9ebc9f3](https://github.com/Danderson123/Amira/commit/9ebc9f31059c589503fb77164c0c171e1627efb4))
* add minimal instructions to README ([570e8af](https://github.com/Danderson123/Amira/commit/570e8afc135d79627c164c7b96e1a33815789b90))
* add pandas as dependency ([5569a2e](https://github.com/Danderson123/Amira/commit/5569a2ec792ae9dacc330a093690516db047409f))
* add sourmash as a dependency ([edc19a0](https://github.com/Danderson123/Amira/commit/edc19a0e75878f113988831a5aa270a126d11552))
* add sourmash as a dependency ([11a4b4d](https://github.com/Danderson123/Amira/commit/11a4b4d43d8a86e41209b40734f796bc89556f7f))
* apply allele filters to genes from short reads ([e3fdc4f](https://github.com/Danderson123/Amira/commit/e3fdc4f116b7550b21aac11623c61a1e2a33619d))
* avoids undercalling genes due to collapsed paths ([02787dd](https://github.com/Danderson123/Amira/commit/02787ddd1955e3cc1a32753724d00063f9d3ce25))
* bug fix for getting the genes from a list of nodes ([6114f75](https://github.com/Danderson123/Amira/commit/6114f75638685881b884bfc6ca8c0396f4b53275))
* bug fix for missing gene due to improper sorting ([3b79f5a](https://github.com/Danderson123/Amira/commit/3b79f5a983c888bc90e1b77259b7527edddb8be5))
* bug in bubble popping that was missing bubbles ([7b1de79](https://github.com/Danderson123/Amira/commit/7b1de79f1a7c3075545d3ea792da315576c0e4f0))
* correct relative import of scripts ([9986fc8](https://github.com/Danderson123/Amira/commit/9986fc86d6455e78ac99f88b4463500ea1684bb1))
* de novo check gene coverages with minimap2 and do not allow racon to trim AMR alleles ([709e0b0](https://github.com/Danderson123/Amira/commit/709e0b0f2896d8d39362a2e2a91c9c6c88ff0564))
* edge case fix for assigning reads to paths and adds a component ID column to the amira tsv ([83cc331](https://github.com/Danderson123/Amira/commit/83cc33198789180751814fedbb71bee42776b6bf))
* edge case in AMR path finding ([95889a8](https://github.com/Danderson123/Amira/commit/95889a8ebb590fcb7b84a40f1dd6a3d38e6dd77d))
* find genes missing from the graph ([108d1df](https://github.com/Danderson123/Amira/commit/108d1dffe62e9dccd0f062e5d59edf906100d01f))
* find genes missing from the graph ([6713e4d](https://github.com/Danderson123/Amira/commit/6713e4df21f7f3e2b83712d4e073319395db651a))
* fix a bunch of bugs that messed with evaluation in real data ([a4e76c9](https://github.com/Danderson123/Amira/commit/a4e76c9a141e4e2cb3ad7410b53e914dc139a2ac))
* fix for known bug where two different AMR paths overlap, tweaks params to filer out AMR genes and prevents k from being increased when the coverage is too low ([70506d1](https://github.com/Danderson123/Amira/commit/70506d13d8f74e085da6859c0889a853ec486024))
* fixes bug where AMR blocks are missed if they start or end at a junction ([c5c0d46](https://github.com/Danderson123/Amira/commit/c5c0d46346e4b53711e523d7251d4bd42771dbd4))
* improper handling of edge cases where duplicate genes occur on a read we are trying to correct and/or in a path that we are trying to correct to. ([6973ecb](https://github.com/Danderson123/Amira/commit/6973ecb492fefe642f1c0699e364f9a64d306168))
* improves modularity of main function and allows JSON inputs ([2ed1347](https://github.com/Danderson123/Amira/commit/2ed1347a678a8fbaf29e7348a5b137208d2f69d7))
* improves modularity of main function and allows JSON inputs ([02d15b0](https://github.com/Danderson123/Amira/commit/02d15b0b2ed2374aaf6711ceb6f04c5f5ddf2552))
* increase coverage threshold for filtering and allow replacement of AMR genes in low coverage paths ([a80d980](https://github.com/Danderson123/Amira/commit/a80d980650c5e2f1dfe3fd65b7d65244f9d2740f))
* limit allele coverage to 100 ([06ae3e2](https://github.com/Danderson123/Amira/commit/06ae3e2364f2eb74b41fc34d2202d4fd45cebc93))
* minor memory optimisation for minimizer extraction ([3dde066](https://github.com/Danderson123/Amira/commit/3dde0668f48565fc28f1582c7b73d3bd292f694b))
* missing AMR path ([8e94b44](https://github.com/Danderson123/Amira/commit/8e94b440d52a6af0ca8af5f6dae77b96a69d0dc9))
* missing AMR path ([b0d64fd](https://github.com/Danderson123/Amira/commit/b0d64fdb48cc95e7a5e43685cf8d1fd95047748d))
* output empty dataframe when no AMR genes are present ([6fef804](https://github.com/Danderson123/Amira/commit/6fef8044c212892522685b8b886f5b9ed1d8f553))
* prevent removal of full components and AMR nodes in dead ends ([deea859](https://github.com/Danderson123/Amira/commit/deea859a3b24c923099d66a4db570f9a1ddc7b17))
* reduce AMR allele coverage threshold to 85% ([92cca4f](https://github.com/Danderson123/Amira/commit/92cca4fee682745f5cdb0f3205e36b99079620be))
* relative import location of fastq functions ([4e74456](https://github.com/Danderson123/Amira/commit/4e7445678472f4bff9d50e06ac5547b6df0e20e5))
* remove legacy script ([3094671](https://github.com/Danderson123/Amira/commit/3094671df13d64763c7dc03cd3bd042b73d8687f))
* remove redundant print ([d9908df](https://github.com/Danderson123/Amira/commit/d9908df5650452562c86d4517a3b3075eff2b6fe))
* remove threshold for initial AMR path finding ([551e090](https://github.com/Danderson123/Amira/commit/551e09096fef4d50abca958651941eb434b08472))
* return empty result when no AMR genes found ([659fceb](https://github.com/Danderson123/Amira/commit/659fcebdc917413b80f938418f97cdc6917e489b))
* sort ref alleles first by length then by similarity ([e056b56](https://github.com/Danderson123/Amira/commit/e056b569699c9bccaaaeff6d71b9ef48db889cf9))
* tweaks parameters to filter AMR genes and ensures subsampling reads does not remove AMR genes completely ([e2ff182](https://github.com/Danderson123/Amira/commit/e2ff1820453fda0feff70a4c26fd70b584012d98))
* undefined variable fix and pre-commit reformatting ([9e9c8e0](https://github.com/Danderson123/Amira/commit/9e9c8e0b59f26d86b088b5dd2a1a366aa9615368))
* update to latest pysam ([beebe6d](https://github.com/Danderson123/Amira/commit/beebe6d75544721884bbff944ed85e1ee03fd591))
* use reads to define start and stop points of paths through AMR blocks ([88dcf51](https://github.com/Danderson123/Amira/commit/88dcf51648986b51c52ae5f427515ed3f7f51c5f))


### Styles

* black reformatting ([c2c2ed0](https://github.com/Danderson123/Amira/commit/c2c2ed0195760ed1d4c5aef8b4eadf303fa80f57))
* black reformatting ([3502cfc](https://github.com/Danderson123/Amira/commit/3502cfcf988f88937e4c7e61acdbb529d00b23b5))
* black refotmatting recommendations ([b41c635](https://github.com/Danderson123/Amira/commit/b41c635b3e4468a74ded92f1d9e05b9fbb1500b2))
* fix all flake8 style recommendations ([13f66f1](https://github.com/Danderson123/Amira/commit/13f66f1e802c8ef60c8c6dfc2d71a302ba38f96a))
* flake8 formatting ([bd28289](https://github.com/Danderson123/Amira/commit/bd282899e17bb4ed58744805425dad9219cd62b2))
* flake8 formatting ([0d477de](https://github.com/Danderson123/Amira/commit/0d477de340049ef9a340c06930d8f94f5c092f93))
* flake8 formatting ([436d5d9](https://github.com/Danderson123/Amira/commit/436d5d9f8a5a3828dfb1efd04f9d92b8b520dbbb))
* flake8 formatting recommendations ([d7d37b8](https://github.com/Danderson123/Amira/commit/d7d37b869eeeae07961b40a72b13a439a457c0c8))
* flake8 formatting recommendations ([b02baec](https://github.com/Danderson123/Amira/commit/b02baec1df8b068f6e641c56560112b803afe763))
* flake8 reformatting ([0f9333c](https://github.com/Danderson123/Amira/commit/0f9333cdd2c6929cfe93b4c41f6e252fceb5d0f3))
* flake8 reformatting ([06a72f2](https://github.com/Danderson123/Amira/commit/06a72f24ff0be0b9fd0c3768b416efeaf7dd8453))
* flake8 reformatting ([c93b288](https://github.com/Danderson123/Amira/commit/c93b28844fa1a34b8f817dcb7e76bdcfb0d45d88))
* flake8 reformatting ([da83389](https://github.com/Danderson123/Amira/commit/da833890d5937b88b8db1c3ccb6bb834a9fda12c))
* isort and black recommended formatting ([51853a6](https://github.com/Danderson123/Amira/commit/51853a6940ed3bd1acc95cd433b9be69d35f502b))
* isort and black recommended formatting ([92efeff](https://github.com/Danderson123/Amira/commit/92efeffa4e3ca5fc9e8d63b5ad2c5afd980f6c77))
* partial typing ([4abdcd4](https://github.com/Danderson123/Amira/commit/4abdcd490e511d6a88f1488bb7fc61911a81e445))
* pre-commit formatting ([c4f18d7](https://github.com/Danderson123/Amira/commit/c4f18d783d18a45c45edf3825a044257c68973f1))
* pre-commit recommended formatting ([8a097d9](https://github.com/Danderson123/Amira/commit/8a097d952e3986e344970535c32a9cdb2de05f03))
* pre-commit reformat ([0d6a1f2](https://github.com/Danderson123/Amira/commit/0d6a1f2fe9301c0fe9d2f6230ac3542480b2bc5b))
* pre-commit reformatting ([e188bb9](https://github.com/Danderson123/Amira/commit/e188bb9791b0414538212a35a89d4745e00fef56))
* pre-commit reformatting ([3ae42d7](https://github.com/Danderson123/Amira/commit/3ae42d767f0575375d0f7343662dd8296ac79fe0))
* reformatting recommendations ([e8f2b57](https://github.com/Danderson123/Amira/commit/e8f2b5737fece61d07235f191eaea2a15f824c73))
* reformatting recommendations ([e523293](https://github.com/Danderson123/Amira/commit/e523293681cf453ebc7730387d986aa15f6cab81))
* reformatting using black ([451a5b4](https://github.com/Danderson123/Amira/commit/451a5b4744099e88f4688fab03f5e1d116130a58))
* remove redundant commented code ([cec5de5](https://github.com/Danderson123/Amira/commit/cec5de5969fd87fdd2ae38df2f6dc29b5d60b06e))
* rename to Amira ([3f4da1c](https://github.com/Danderson123/Amira/commit/3f4da1c1b44888b3a976d4baebc1eb4237aa578c))
* rename to Amira ([7a4a6f5](https://github.com/Danderson123/Amira/commit/7a4a6f58232c2b6c466fa0dff48bd82e109a6a02))


### Tests

* add additional test files ([f606121](https://github.com/Danderson123/Amira/commit/f606121121232e91215b97a66bb392b5d6d80d5c))
* add more complexity in tests ([00ef219](https://github.com/Danderson123/Amira/commit/00ef21908f5f3d614afb87da290fa09af9c85e50))
* add test files ([f3f9a51](https://github.com/Danderson123/Amira/commit/f3f9a51f303ac02e4155068f8dea6c86184bf4b9))
* add unittest for minimizer comparison ([4c4934d](https://github.com/Danderson123/Amira/commit/4c4934d19fc3c35b9fe1af5db5581cb76e9b15bf))
* add unittest for overlapping AMR blocks ([2f69561](https://github.com/Danderson123/Amira/commit/2f69561b3a92a64d314199606fbd7d3501a65430))
* remove redundant test ([cc79dc0](https://github.com/Danderson123/Amira/commit/cc79dc0e8c04ad52058e56a4bb0ec6d3c929a8a6))
* remove test graph ([d308ee3](https://github.com/Danderson123/Amira/commit/d308ee3b75211a964ff24ce12fec1bd932cddd83))
