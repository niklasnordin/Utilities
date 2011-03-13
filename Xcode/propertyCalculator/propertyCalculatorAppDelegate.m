//
//  propertyCalculatorAppDelegate.m
//  propertyCalculator
//
//  Created by Niklas Nordin on 2010-12-30.
//  Copyright 2010 nequam. All rights reserved.
//

#import "propertyCalculatorAppDelegate.h"
#import "db.h"

@implementation propertyCalculatorAppDelegate

@synthesize window;

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
	// Insert code here to initialize your application 
	db *database = [[db alloc] init];
	[database create];
}

@end
