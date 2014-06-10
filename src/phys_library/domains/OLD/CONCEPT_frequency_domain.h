public:

  static int                  get_size();
  static std::vector<double>& get_elements();
  
  static int phase(double frequency);

  static void to_JSN(std::stringstream& ss);

private:

  static std::vector<double>& make_elements();
