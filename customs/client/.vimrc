scriptencoding utf-8
set nocompatible               " be iMproved
filetype off                   " required!

set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()

Plugin 'VundleVim/Vundle.vim'

" My Bundles here:
" original repos on github
Plugin 'Valloric/YouCompleteMe'
Plugin 'scrooloose/syntastic'
Plugin 'scrooloose/nerdcommenter'
Plugin 'SirVer/ultisnips'
Plugin 'vim-jp/cpp-vim'
Plugin 'rhysd/vim-clang-format'
Plugin 'realincubus/vim-clang-refactor'
Plugin 'octol/vim-cpp-enhanced-highlight'
Plugin 'ervandew/supertab'
Plugin 'tpope/vim-fugitive'
Plugin 'tpope/vim-surround'
Plugin 'flazz/vim-colorschemes'
Plugin 'dikiaap/minimalist'
Plugin 'vim-airline/vim-airline'

" vimcountry vimrc
Plugin 'tpope/vim-rails'
Plugin 'tpope/vim-haml'
Plugin 'tpope/vim-endwise'
Plugin 'tpope/vim-ragtag'
Plugin 'tpope/vim-markdown'
Plugin 'tpope/vim-unimpaired'
Plugin 'scrooloose/nerdtree'
Plugin 'tomtom/tlib_vim'
Plugin 'tomtom/tcomment_vim'
Plugin 'tomtom/tselectbuffer_vim'
Plugin 'trapd00r/x11colors.vim'
Plugin 'lilydjwg/colorizer'
Plugin 'shemerey/vim-project'
Plugin 'Twinside/vim-codeoverview'
Plugin 'msanders/snipmate.vim'
Plugin 'Townk/vim-autoclose'

" LaTeX Plugins
Plugin 'vim-latex/vim-latex'

call vundle#end()

" General settings {{{
	filetype on
	filetype plugin indent on
	syntax on

	set encoding=utf-8
	set fileencoding=utf-8

	set title
	set mouse=a
	set hlsearch

	set shortmess=at      " shorten error messages

	set nrformats+=alpha  " in-/decrease letters with C-a/C-x

	set modeline          " enable modelines
	set modelines=5

	set number            " enable line numbers
	set ruler             " enable something
	set cursorline        " enable hiliting of cursor line

	set backspace=2       " backspace over EOL etc.

	set background=dark   " i prefer dark backgrounds

	set hidden            " buffer switching should be quick
	set confirm           " ask instead of just print errors
	set equalalways       " make splits equal size

	set lazyredraw        " don't redraw while executing macros

	set noshowmode        " don't display mode, it's already in the status line

	let mapleader=","
	let maplocalleader=","

	" airline {{{
		let g:airline_powerline_fonts = 1
		let g:airline#extensions#tabline#enabled = 1
	" }}}
	" ycm clang {{{
		let g:clang_c_options = '-std=gnu11'
		let g:clang_cpp_options = '-std=c++11 -stdlib=libc++'
		let g:ycm_confirm_extra_conf = 0 

		let g:ycm_key_list_select_completion = ['<C-n>', '<Down>']
		let g:ycm_key_list_previous_completion = ['<C-p>', '<Up>']
		let g:ycm_server_log_level = "debug"
		let g:SuperTabDefaultCompletionType = '<C-n>'
	" }}}
	" ultisnips {{{
	" better key bindings for UltiSnipsExpandTrigger
		let g:UltiSnipsExpandTrigger = "<tab>"
		let g:UltiSnipsJumpForwardTrigger = "<tab>"
		let g:UltiSnipsJumpBackwardTrigger = "<s-tab>"
		let g:UltiSnipsEditSplit="horizontal"
	" }}}
	" syntastic {{{
		let g:syntastic_always_populate_loc_list = 1
		let g:syntastic_auto_loc_list = 1
		let g:syntastic_check_on_open = 1
		let g:syntastic_check_on_wq = 0
	" }}}
" }}}
" vim-latex settings {{{
	" IMPORTANT: win32 users will need to have 'shellslash' set so that latex
	" can be called correctly.
	set shellslash
	" IMPORTANT: grep will sometimes skip displaying the file name if you
	" search in a singe file. This will confuse Latex-Suite. Set your grep
	" program to always generate a file-name.
	set grepprg=grep\ -nH\ $*
	" The following changes the default filetype back to 'tex':
	let g:tex_flavor='latex'
" }}}
" General Keybinds {{{

  " Set MapLeader
  let mapleader = ","

  " Delete previous word with C-BS
  imap <C-BS> <C-W>

  " Toggle Buffer Selection and list Tag Lists
  map <F2> <Esc>:TSelectBuffer<CR>
  map <F4> <Esc>:TlistToggle<CR>

   
  " Set text wrapping toggles
  nmap <silent> <leader>w :set invwrap<CR>:set wrap?<CR> 
   
  " Set up retabbing on a source file
  nmap <silent> <leader>rr :1,$retab<CR> 
   
  " cd to the directory containing the file in the buffer
  nmap <silent> <leader>cd :lcd %:h<CR> 
   
  " Make the directory that contains the file in the current buffer.
  " This is useful when you edit a file in a directory that doesn't
  " (yet) exist
  nmap <silent> <leader>md :!mkdir -p %:p:h<CR>

  " Increase @revision # by 1
  nmap <silent> <leader>incr /@updated
	" wwwd$"=strftime("%a %d %b %Y")
	" p/@revision
	" $

" }}}
" {{{ Window movement
  nmap <M-h> :winc h<CR>
  nmap <M-j> :winc j<CR>
  nmap <M-k> :winc k<CR>
  nmap <M-l> :winc l<CR>
" }}}
" GUI or no GUI, that's the question {{{
  if has('gui_running')
    set guicursor+=a:blinkon0       " Cursor doesn't blink - it's annoying
    set guioptions-=m               " No Menubar
    set guioptions-=T               " No Toolbar
    set guioptions-=l               " No Scrollbar left
    set guioptions-=L               " No Scrollbar left when split
    set guioptions-=r               " No Scrollbar right
    set guioptions-=r               " No Scrollbar right when split

    set laststatus=2                " always show statusline

    " set gfn=Pragmata\ 6.5
    set gfn=Neep\ Medium\ Semi-Condensed\ 9
    " set gfn=Mensch\ 7

    set lines=40                    " Height
    set columns=85                  " Width

    colorscheme minimalist

  else
    colorscheme minimalist
  endif
" }}}
" Status line {{{
  set laststatus=2      " always show statusline

  " Generic Statusline {{{
  function! SetStatus()
    setl statusline+=
          \%1*\ %f
          \%H%M%R%W%7*\ ┃
          \%2*\ %Y\ <<<\ %{&ff}%7*\ ┃
  endfunction

  function! SetRightStatus()
    setl statusline+=
          \%5*\ %{StatusFileencoding()}%7*\ ┃
          \%5*\ %{StatusBuffersize()}%7*\ ┃
          \%=%<%7*\ ┃
          \%5*\ %{StatusWrapON()}
          \%6*%{StatusWrapOFF()}\ %7*┃
          \%5*\ %{StatusInvisiblesON()}
          \%6*%{StatusInvisiblesOFF()}\ %7*┃
          \%5*\ %{StatusExpandtabON()}
          \%6*%{StatusExpandtabOFF()}\ %7*┃
          \%5*\ w%{StatusTabstop()}\ %7*┃
          \%3*\ %l,%c\ >>>\ %P
          \\ 
  endfunction " }}}

  " Update when leaving Buffer {{{
  function! SetStatusLeaveBuffer()
    setl statusline=""
    call SetStatus()
  endfunction
  au BufLeave * call SetStatusLeaveBuffer() " }}}

  " Update when switching mode {{{
  function! SetStatusInsertMode(mode)
    setl statusline=%4*
    if a:mode == 'i'
      setl statusline+=\ Einfügen\ ◥
    elseif a:mode == 'r'
      setl statusline+=\ Ersetzen\ ◥
    elseif a:mode == 'normal'
      setl statusline+=\ \ ◥
    endif
    call SetStatus()
    call SetRightStatus()
  endfunction

  au VimEnter     * call SetStatusInsertMode('normal')
  au InsertEnter  * call SetStatusInsertMode(v:insertmode)
  au InsertLeave  * call SetStatusInsertMode('normal')
  au BufEnter     * call SetStatusInsertMode('normal') " }}}

  " Some Functions shamelessly ripped and modified from Cream
  "fileencoding (three characters only) {{{
  function! StatusFileencoding()
    if &fileencoding == ""
      if &encoding != ""
        return &encoding
      else
        return " -- "
      endif
    else
      return &fileencoding
    endif
  endfunc " }}}
  " &expandtab {{{
  function! StatusExpandtabON()
    if &expandtab == 0
      return "tabs"
    else
      return ""
    endif
  endfunction "
  function! StatusExpandtabOFF()
    if &expandtab == 0
      return ""
    else
      return "tabs"
    endif
  endfunction " }}}
  " tabstop and softtabstop {{{
  function! StatusTabstop()

    " show by Vim option, not Cream global (modelines)
    let str = "" . &tabstop
    " show softtabstop or shiftwidth if not equal tabstop
    if   (&softtabstop && (&softtabstop != &tabstop))
    \ || (&shiftwidth  && (&shiftwidth  != &tabstop))
      if &softtabstop
        let str = str . ":sts" . &softtabstop
      endif
      if &shiftwidth != &tabstop
        let str = str . ":sw" . &shiftwidth
      endif
    endif
    return str

  endfunction " }}}
  " Buffer Size {{{
  function! StatusBuffersize()
    let bufsize = line2byte(line("$") + 1) - 1
    " prevent negative numbers (non-existant buffers)
    if bufsize < 0
      let bufsize = 0
    endif
    " add commas
    let remain = bufsize
    let bufsize = ""
    while strlen(remain) > 3
      let bufsize = "," . strpart(remain, strlen(remain) - 3) . bufsize
      let remain = strpart(remain, 0, strlen(remain) - 3)
    endwhile
    let bufsize = remain . bufsize
    " too bad we can't use "¿" (nr2char(1068)) :)
    let char = "b"
    return bufsize . char
  endfunction " }}}
  " Show Invisibles {{{
  function! StatusInvisiblesON()
    "if exists("g:LIST") && g:LIST == 1
    if &list
      if     &encoding == "latin1"
        return "¶"
      elseif &encoding == "utf-8"
        return "¶"
      else
        return "$"
      endif
    else
      return ""
    endif
  endfunction
  function! StatusInvisiblesOFF()
    "if exists("g:LIST") && g:LIST == 1
    if &list
      return ""
    else
      if     &encoding == "latin1"
        return "¶"
      elseif &encoding == "utf-8"
        return "¶"
      else
        return "$"
      endif
    endif
  endfunction " }}}
  " Wrap Enabled {{{
  function! StatusWrapON()
    if &wrap
      return "wrap"
    else
      return ""
    endif
  endfunction
  function! StatusWrapOFF()
    if &wrap
      return ""
    else
      return "wrap"
    endif
  endfunction
  " }}}
" }}}
" Tabstops {{{
  set tabstop=2
  set shiftwidth=2
  set softtabstop=2
  set autoindent
  set smartindent
  set expandtab
" }}}
" Invisibles {{{
"  set listchars=tab:>\ ,eol:<
"  set list
"  nmap <silent> <F5> :set list!<CR>
" }}}
" Tabstops {{{
  set tabstop=2
  set shiftwidth=2
  set softtabstop=2
  set autoindent
  set smartindent
  set expandtab
" }}}
" Invisibles {{{
"  set listchars=tab:>\ ,eol:<
"  set list
"  nmap <silent> <F5> :set list!<CR>
" }}}
" Folds {{{
  set foldmethod=marker
  set foldcolumn=1
  " au BufWinLeave * mkview
  " au BufWinEnter * silent loadview
" }}}
" Pairings {{{
  set showmatch
" }}}
" Margins {{{
  set scrolloff=5
  set sidescroll=5
" }}}
" Search {{{
  set incsearch
  set ignorecase

  " Toggle that stupid highlight search
  nmap <silent> ,n :set invhls<CR>:set hls?<CR> 
" }}}
" Backup files {{{
  set nobackup
  set nowb
  set noswapfile
" }}}
" Completion {{{
  set wildmenu
  set wildmode=longest,full,list

  set ofu=syntaxcomplete#Complete
" }}}
" NERDTree {{{
  map <F3> :NERDTreeToggle<CR>

  let NERDTreeChDirMode = 2
  let NERDTreeShowBookmarks = 1
" }}}
" Wrapping {{{
  set linebreak
  set showbreak=↳\ 
" toggle wrapping
  nmap <silent> <F12> :let &wrap = !&wrap<CR>
" }}}
" RagTag {{{
  imap <M-O> <Esc>o
  imap <C-J> <Down>
  let g:ragtag_global_maps = 1

  imap <C-Space> <C-X><Space>
  imap <C-CR> <C-X><CR>
" }}}
" 'Bubbling' {{{
  nmap <C-up> [e
  nmap <C-down> ]e
  vmap <C-up> [egv
  vmap <C-down> ]egv
" }}}
" Formatting with Par (gqip) {{{
  set formatprg=par\ -req
  nmap <F9> gqip
" }}}
" Pasting {{{
  set paste
  nnoremap p ]p
  nnoremap <c-p> p
" }}}
" Macros {{{
  " Execute macro "q" with space
  nmap <Space> @q
  " Map @ to + for more comfortable macros on DE kb layout
  nmap + @
" }}}