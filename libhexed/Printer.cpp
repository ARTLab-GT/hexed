#include <hexed/Printer.hpp>
#include <hexed/utils.hpp>

namespace hexed {

Compound_printer::Compound_printer(std::vector<std::shared_ptr<Printer>> p) : printers{p} {}

void Compound_printer::operator()(std::string message, bool emph, bool replace) {
  for (auto& printer : printers) (*printer)(message, emph, replace);
}

const Stream_printer::Format Stream_printer::default_format;

Stream_printer::Stream_printer(std::ostream& stream, Format emph_format)
: _stream{stream}, _format_code(""), _reset_code("")
{
  if (!(emph_format.type == unspecified_type && emph_format.color == unspecified_color)) {
    std::string type_format = emph_format.type == unspecified_type ? "" : std::to_string(emph_format.type);
    if (emph_format.color == unspecified_color) emph_format.color = default_color;
    _format_code = format_str(100, "\x1b[%s;%i%im", type_format.c_str(), 3 + 6*emph_format.light + emph_format.background, emph_format.color);
    _reset_code = "\x1b[0m";
  }
}

void Stream_printer::operator()(std::string message, bool emph, bool replace) {
  #pragma omp critical
  {
    if (replace) _stream << "\x1b[G\x1b[K";
    _stream << (emph ? _format_code : "") << message << (emph ? _reset_code : "") << std::flush;
  }
}

Compound_printer printers::info(std::vector<std::shared_ptr<Printer>> {
  std::make_shared<Stream_printer>(std::cout, Stream_printer::Format {.type = Stream_printer::bold,
                                                                      .color = Stream_printer::green})
});

Compound_printer printers::warn(std::vector<std::shared_ptr<Printer>> {
  std::make_shared<Stream_printer>(std::cerr, Stream_printer::Format {.color = Stream_printer::yellow, .light = true})
});

Compound_printer printers::error(std::vector<std::shared_ptr<Printer>> {
  std::make_shared<Stream_printer>(std::cerr, Stream_printer::Format {.type = Stream_printer::bold,
                                                                      .color = Stream_printer::red})
});

Printer_set::Printer_set() {
  info.printers.emplace_back(std::make_shared<Stream_printer>(std::cout, Stream_printer::Format{.type = Stream_printer::bold,
                                                                                                .color = Stream_printer::green}));
  warn.printers.emplace_back(std::make_shared<Stream_printer>(std::cerr, Stream_printer::Format{.color = Stream_printer::yellow, .light = true}));
  error.printers.emplace_back(std::make_shared<Stream_printer>(std::cerr, Stream_printer::Format{.type = Stream_printer::bold,
                                                                                                 .color = Stream_printer::red}));
}

}
