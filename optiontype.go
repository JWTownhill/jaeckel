package jaeckel

// OptionType represents the type of option (call or put).
type OptionType int

const (
	// Call option.
	Call OptionType = 1
	// Put option.
	Put OptionType = -1
)
