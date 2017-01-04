scriptencoding utf-8

set nocompatible               " be iMproved
filetype off                   " required!

set rtp+=/opt/vim/bundle/vundle
call vundle#rc("/opt/vim/bundle/")

Bundle 'gmarik/vundle'

" My Bundles here:
" original repos on github
Bundle 'Valloric/YouCompleteMe'
Bundle 'scrooloose/syntastic'
Bundle 'SirVer/ultisnips'
Bundle 'vim-jp/cpp-vim'
Bundle 'rhysd/vim-clang-format'
Bundle 'realincubus/vim-snippets'
Bundle 'realincubus/vim-clang-refactor'
Bundle 'octol/vim-cpp-enhanced-highlight'
Bundle 'ervandew/supertab'


let g:clang_c_options = '-std=gnu11'
let g:clang_cpp_options = '-std=c++11 -stdlib=libc++'
let g:ycm_confirm_extra_conf = 0 

filetype plugin indent on     " required!

set shiftwidth=2
set softtabstop=2
set tabstop=2
set number
set autoindent

let g:ycm_key_list_select_completion = ['<C-n>', '<Down>']
let g:ycm_key_list_previous_completion = ['<C-p>', '<Up>']
let g:ycm_server_log_level = "debug"
let g:SuperTabDefaultCompletionType = '<C-n>'

" better key bindings for UltiSnipsExpandTrigger
let g:UltiSnipsExpandTrigger = "<tab>"
let g:UltiSnipsJumpForwardTrigger = "<tab>"
let g:UltiSnipsJumpBackwardTrigger = "<s-tab>"

let g:UltiSnipsEditSplit="horizontal"

syntax on

set nobackup
set nowritebackup

colorscheme badwolf

set mouse-=a
set background=light


